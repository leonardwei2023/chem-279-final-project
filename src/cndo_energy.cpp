/**
 * cndo_energy.cpp
 *
 * CNDO/2 SCF total energy calculator.
 * Called by CNDOEngine via shell command for finite-difference Hessian.
 *
 * Theory: Lectures 13-14 (CHEM 279, Mayank Agrawal)
 *
 * CNDO/2 Fock matrix (Lec 13, slide 7):
 *
 *   Diagonal (orbital mu on atom A):
 *     F_mumu = -½(I_mu + A_mu)
 *              + [P_AA^tot - Z_A - (P_mumu^alpha - ½)] * gamma_AA
 *              + sum_{B≠A} [P_BB^tot - Z_B] * gamma_AB
 *
 *   Off-diagonal (mu on A, nu on B, A≠B):
 *     F_munu = -½(beta_A + beta_B) * S_munu - P_munu^alpha * gamma_AB
 *
 * For closed-shell RHF: P^alpha_munu = P^total_munu / 2
 *   where P^total = 2 * C_occ * C_occ^T
 *
 * Density matrix (Lec 13, slide 8):
 *   P^alpha_munu = sum_{i=1}^{n_occ} C_mui * C_nui     (n_occ alpha electrons)
 *   P^tot_munu   = 2 * P^alpha_munu
 *
 * Total energy (Lec 12, standard CNDO/2):
 *   E = ½ sum_{mu,nu} P^tot_munu * (H^core_munu + F_munu)
 *     + sum_{A<B} Z_A * Z_B * gamma_AB
 *
 * Core Hamiltonian (Beveridge & Pople CNDO/2):
 *   H^core_mumu = -½(I_mu + A_mu) + (Z_A - ½)*gamma_AA - sum_{B≠A} Z_B*gamma_AB
 *   H^core_munu = -½(beta_A + beta_B)*S_munu   (A≠B)
 *
 * Gamma integrals (Lec 14, slide 4 - Mataga-Nishimoto approximation):
 *   gamma_AA  = I_s - A_s   (on-site, from ionization potential / electron affinity)
 *   gamma_AB  = 1 / (R_AB + 2/(gamma_AA + gamma_BB))   (two-center)
 *
 * Overlap integrals: STO-3G contracted Gaussians (Lec 7-8, Lec 10)
 *
 * Valence basis (CNDO/2 is valence-only, Lec 13 slide 9):
 *   H:           1s  (1 orbital)
 *   C,N,O,F,Cl:  2s, 2px, 2py, 2pz  (4 orbitals — NO core 1s)
 *
 * Atomic parameters (Beveridge & Pople / Pople & Segal, CNDO/2):
 *   H:  I_s=13.06 eV, A_s=3.69 eV,  beta=-9.0 eV,  zeta=1.2
 *   C:  I_s=14.051, A_s=5.572, I_p=5.572, A_p=2.275, beta=-21.0, zeta=1.625
 *   N:  I_s=19.316, A_s=7.275, I_p=7.275, A_p=2.862, beta=-25.0, zeta=1.950
 *   O:  I_s=25.390, A_s=9.111, I_p=9.111, A_p=3.524, beta=-31.0, zeta=2.275
 *   F:  I_s=32.272, A_s=11.08, I_p=11.08, A_p=4.270, beta=-39.0, zeta=2.425
 *   Cl: I_s=24.350, A_s=3.890, I_p=3.890, A_p=1.430, beta=-15.0, zeta_s=3.500, zeta_p=2.033
 */

#include "molecule.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ─── Constants ───────────────────────────────────────────────────────────────
static constexpr double EV_TO_HA    = 0.03674930;     // 1 eV = 0.036749 Hartree
static constexpr double ANG_TO_BOHR = 1.8897261245;   // 1 Angstrom = 1.8897 Bohr
static constexpr double PI          = 3.14159265358979323846;

// ─── STO-3G Gaussian primitives ──────────────────────────────────────────────
// Source: Hehre, Stewart, Pople (1969); referenced in Lec 10
// For orbital with STO exponent zeta, actual Gaussian exponent = alpha_j * zeta^2
struct GTO { double alpha, coeff; };

static const GTO STO3G_1s[3] = {
    {2.2277, 0.1543}, {0.4058, 0.5353}, {0.1098, 0.4446}
};
static const GTO STO3G_2s[3] = {
    {2.5820, -0.0999}, {0.1567, 0.3995}, {0.0611, 0.7002}
};
static const GTO STO3G_2p[3] = {
    {2.5820, 0.1559}, {0.1567, 0.6077}, {0.0611, 0.3920}
};

// ─── Atomic parameters ───────────────────────────────────────────────────────
struct AtomParams {
    int    z_valence;    // valence electron count (CNDO/2 treats valence only)
    int    n_orbitals;   // 1 for H, 4 for 2nd-row elements
    double I_s, A_s;     // s ionization potential and electron affinity (eV)
    double I_p, A_p;     // p ionization potential and electron affinity (eV)
    double beta;         // resonance parameter beta (eV)
    double zeta_s;       // STO exponent for s orbital
    double zeta_p;       // STO exponent for p orbital
};

static AtomParams get_params(const std::string& sym) {
    //                    Zv norb  I_s     A_s    I_p    A_p    beta   zeta_s zeta_p
    if (sym == "H")  return {1, 1, 13.06,  3.69,  0.0,   0.0,   -9.0,  1.200, 0.000};
    if (sym == "C")  return {4, 4, 14.051, 5.572, 5.572, 2.275, -21.0, 1.625, 1.625};
    if (sym == "N")  return {5, 4, 19.316, 7.275, 7.275, 2.862, -25.0, 1.950, 1.950};
    if (sym == "O")  return {6, 4, 25.390, 9.111, 9.111, 3.524, -31.0, 2.275, 2.275};
    if (sym == "F")  return {7, 4, 32.272, 11.08, 11.08, 4.270, -39.0, 2.425, 2.425};
    if (sym == "Cl") return {7, 4, 24.350, 3.890, 3.890, 1.430, -15.0, 3.500, 2.033};
    throw std::runtime_error("cndo_energy: unsupported element: " + sym);
}

// ─── Orbital type ────────────────────────────────────────────────────────────
enum OrbType { ORB_1S, ORB_2S, ORB_2PX, ORB_2PY, ORB_2PZ };

struct OrbInfo {
    int     atom_idx;
    OrbType type;
    double  zeta;
    double  half_IP;   // ½(I_mu + A_mu) in Hartree
};

// Build valence orbital list — CNDO/2 excludes core 1s for heavy atoms
static std::vector<OrbInfo> make_orbitals(
    const Molecule& mol,
    const std::vector<AtomParams>& params)
{
    std::vector<OrbInfo> orbs;
    auto syms = mol.get_symbols();
    int n = mol.get_num_atoms();

    for (int A = 0; A < n; A++) {
        const AtomParams& p = params[A];

        if (p.n_orbitals == 1) {
            // Hydrogen: one valence orbital (1s)
            double half_s = 0.5 * (p.I_s + p.A_s) * EV_TO_HA;
            orbs.push_back({A, ORB_1S, p.zeta_s, half_s});
        } else {
            // 2nd-row elements: valence basis = 2s, 2px, 2py, 2pz
            // Core 1s is NOT included in CNDO/2 (Lec 13, slide 9)
            double half_s = 0.5 * (p.I_s + p.A_s) * EV_TO_HA;
            double half_p = 0.5 * (p.I_p + p.A_p) * EV_TO_HA;
            orbs.push_back({A, ORB_2S,  p.zeta_s, half_s});
            orbs.push_back({A, ORB_2PX, p.zeta_p, half_p});
            orbs.push_back({A, ORB_2PY, p.zeta_p, half_p});
            orbs.push_back({A, ORB_2PZ, p.zeta_p, half_p});
        }
    }
    return orbs;
}

// ─── STO-3G overlap integrals ─────────────────────────────────────────────────
// s-s overlap between two primitive Gaussians centered at Ra, Rb
static double prim_ss(double a, double b,
                      const double Ra[3], const double Rb[3])
{
    double dx = Ra[0]-Rb[0], dy = Ra[1]-Rb[1], dz = Ra[2]-Rb[2];
    double R2 = dx*dx + dy*dy + dz*dz;
    return std::pow(PI/(a+b), 1.5) * std::exp(-a*b/(a+b)*R2);
}

// p_i (component pi) on center Ra  vs  s on center Rb
static double prim_ps(int pi, double a, double b,
                      const double Ra[3], const double Rb[3])
{
    double ab = a + b;
    double P[3];
    for (int k = 0; k < 3; k++) P[k] = (a*Ra[k] + b*Rb[k]) / ab;
    double dR2 = 0.0;
    for (int k = 0; k < 3; k++) { double d=Ra[k]-Rb[k]; dR2+=d*d; }
    double PmA = P[pi] - Ra[pi];
    return PmA * std::pow(PI/ab, 1.5) * std::exp(-a*b/ab*dR2);
}

// p_i on Ra  vs  p_j on Rb
static double prim_pp(int pi, int pj, double a, double b,
                      const double Ra[3], const double Rb[3])
{
    double ab = a + b;
    double P[3];
    for (int k = 0; k < 3; k++) P[k] = (a*Ra[k] + b*Rb[k]) / ab;
    double dR2 = 0.0;
    for (int k = 0; k < 3; k++) { double d=Ra[k]-Rb[k]; dR2+=d*d; }
    double PmA = P[pi] - Ra[pi];
    double PmB = P[pj] - Rb[pj];
    double s   = std::pow(PI/ab, 1.5) * std::exp(-a*b/ab*dR2);
    return (PmA*PmB + (pi==pj ? 0.5/ab : 0.0)) * s;
}

// Contracted STO-3G overlap integral between two orbitals
static double sto3g_overlap(const OrbInfo& mu, const double* Rmu,
                             const OrbInfo& nu, const double* Rnu)
{
    auto gtos = [](OrbType t) -> const GTO* {
        if (t == ORB_1S)  return STO3G_1s;
        if (t == ORB_2S)  return STO3G_2s;
        return STO3G_2p;  // 2PX, 2PY, 2PZ
    };

    const GTO* ga = gtos(mu.type);
    const GTO* gb = gtos(nu.type);

    bool is_p_a = (mu.type==ORB_2PX||mu.type==ORB_2PY||mu.type==ORB_2PZ);
    bool is_p_b = (nu.type==ORB_2PX||nu.type==ORB_2PY||nu.type==ORB_2PZ);
    int  pi_a   = (mu.type==ORB_2PX)?0:(mu.type==ORB_2PY)?1:2;
    int  pi_b   = (nu.type==ORB_2PX)?0:(nu.type==ORB_2PY)?1:2;

    double S = 0.0;
    for (int j = 0; j < 3; j++) {
        double alpha_a = ga[j].alpha * mu.zeta * mu.zeta;
        for (int k = 0; k < 3; k++) {
            double alpha_b = gb[k].alpha * nu.zeta * nu.zeta;
            double s_prim;
            if      (!is_p_a && !is_p_b) s_prim = prim_ss(alpha_a, alpha_b, Rmu, Rnu);
            else if ( is_p_a && !is_p_b) s_prim = prim_ps(pi_a, alpha_a, alpha_b, Rmu, Rnu);
            else if (!is_p_a &&  is_p_b) s_prim = prim_ps(pi_b, alpha_b, alpha_a, Rnu, Rmu);
            else                          s_prim = prim_pp(pi_a, pi_b, alpha_a, alpha_b, Rmu, Rnu);
            S += ga[j].coeff * gb[k].coeff * s_prim;
        }
    }
    return S;
}

// ─── Gamma integrals ─────────────────────────────────────────────────────────
// On-site: gamma_AA = I_s - A_s  (Lec 14 slide 4, Mataga-Nishimoto)
static double gamma_onsite(const AtomParams& p) {
    return (p.I_s - p.A_s) * EV_TO_HA;
}

// Two-center: gamma_AB = 1 / (R_AB + 2/(gamma_AA + gamma_BB))
static double gamma_two_center(double R_bohr, double gAA, double gBB) {
    return 1.0 / (R_bohr + 2.0 / (gAA + gBB));
}

// ─── CNDO/2 SCF total energy ──────────────────────────────────────────────────
double cndo2_total_energy(const Molecule& mol) {
    int n_atoms = mol.get_num_atoms();
    auto syms   = mol.get_symbols();
    auto coords = mol.get_coordinates();   // Angstrom

    // Atomic parameters and electron count
    std::vector<AtomParams> params(n_atoms);
    int n_elec = 0;
    for (int A = 0; A < n_atoms; A++) {
        params[A] = get_params(syms[A]);
        n_elec   += params[A].z_valence;
    }
    if (n_elec % 2 != 0) {
        throw std::runtime_error(
            "cndo_energy: odd electron count — open shell not supported.");
    }
    int n_occ = n_elec / 2;   // number of doubly-occupied orbitals

    // Atomic positions in Bohr
    std::vector<std::array<double,3>> R(n_atoms);
    for (int A = 0; A < n_atoms; A++) {
        R[A][0] = coords[3*A]   * ANG_TO_BOHR;
        R[A][1] = coords[3*A+1] * ANG_TO_BOHR;
        R[A][2] = coords[3*A+2] * ANG_TO_BOHR;
    }

    // Valence orbital list (H: 1 orbital, heavy atoms: 4 orbitals)
    auto orbs  = make_orbitals(mol, params);
    int  n_orb = static_cast<int>(orbs.size());

    // Beta parameters per orbital (resonance integral, Hartree)
    std::vector<double> beta_orb(n_orb);
    for (int mu = 0; mu < n_orb; mu++) {
        beta_orb[mu] = params[orbs[mu].atom_idx].beta * EV_TO_HA;
    }

    // On-site and two-center gamma integrals
    std::vector<double> gAA(n_atoms);
    for (int A = 0; A < n_atoms; A++) gAA[A] = gamma_onsite(params[A]);

    Eigen::MatrixXd gamma(n_atoms, n_atoms);
    for (int A = 0; A < n_atoms; A++) {
        gamma(A, A) = gAA[A];
        for (int B = A+1; B < n_atoms; B++) {
            double dx = R[A][0]-R[B][0], dy = R[A][1]-R[B][1], dz = R[A][2]-R[B][2];
            double Rab = std::sqrt(dx*dx + dy*dy + dz*dz);
            double g   = gamma_two_center(Rab, gAA[A], gAA[B]);
            gamma(A,B) = g; gamma(B,A) = g;
        }
    }

    // Overlap matrix S (STO-3G contracted Gaussians)
    Eigen::MatrixXd S = Eigen::MatrixXd::Identity(n_orb, n_orb);
    for (int mu = 0; mu < n_orb; mu++) {
        int A = orbs[mu].atom_idx;
        for (int nu = mu+1; nu < n_orb; nu++) {
            int B = orbs[nu].atom_idx;
            double s = sto3g_overlap(orbs[mu], R[A].data(),
                                     orbs[nu], R[B].data());
            S(mu,nu) = s; S(nu,mu) = s;
        }
    }

    // Core Hamiltonian H^core (density-independent part of F)
    // H^core_mumu = -½(I+A) + (Z_A-½)*gamma_AA - sum_{B!=A} Z_B*gamma_AB
    // H^core_munu = -½(beta_A+beta_B)*S_munu   (A!=B, ZDO: same-atom off-diag = 0)
    Eigen::MatrixXd H_core = Eigen::MatrixXd::Zero(n_orb, n_orb);
    for (int mu = 0; mu < n_orb; mu++) {
        int A = orbs[mu].atom_idx;
        double h = -orbs[mu].half_IP;
        h += (params[A].z_valence - 0.5) * gamma(A, A);
        for (int B = 0; B < n_atoms; B++) {
            if (B != A) h -= params[B].z_valence * gamma(A, B);
        }
        H_core(mu, mu) = h;

        for (int nu = mu+1; nu < n_orb; nu++) {
            int B = orbs[nu].atom_idx;
            if (A != B) {
                double hoff = -0.5*(beta_orb[mu]+beta_orb[nu])*S(mu,nu);
                H_core(mu,nu) = hoff; H_core(nu,mu) = hoff;
            }
            // same atom: ZDO → 0 (already zero from Eigen::Zero)
        }
    }

    // Nuclear repulsion (CNDO/2 ZDO form):
    // E_nuc = sum_{A<B} Z_A * Z_B * gamma_AB
    double E_nuc = 0.0;
    for (int A = 0; A < n_atoms; A++)
        for (int B = A+1; B < n_atoms; B++)
            E_nuc += params[A].z_valence * params[B].z_valence * gamma(A,B);

    // ── SCF iterations ────────────────────────────────────────────────────────
    // Initialize: P^total = 0 (standard starting guess)
    Eigen::MatrixXd P_tot = Eigen::MatrixXd::Zero(n_orb, n_orb);

    const int    MAX_ITER  = 300;
    const double CONV_DENS = 1.0e-8;   // RMS change in density matrix
    const double CONV_ENER = 1.0e-10;  // change in total energy

    double E_old = 1.0e30;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        // Per-atom total electron populations: P_AA^tot = sum_{mu on A} P^tot_mumu
        std::vector<double> P_AA(n_atoms, 0.0);
        for (int mu = 0; mu < n_orb; mu++)
            P_AA[orbs[mu].atom_idx] += P_tot(mu, mu);

        // Alpha-spin density: P^alpha = P^total / 2  (closed shell)
        Eigen::MatrixXd P_alpha = 0.5 * P_tot;

        // Build Fock matrix (Lec 13, slide 7)
        // F_mumu = -½(I+A) + [P_AA^tot - Z_A - (P^alpha_mumu - ½)]*gamma_AA
        //        + sum_{B!=A} [P_BB^tot - Z_B]*gamma_AB
        // F_munu = -½(beta_A+beta_B)*S_munu - P^alpha_munu*gamma_AB  (A!=B)
        Eigen::MatrixXd F(n_orb, n_orb);

        for (int mu = 0; mu < n_orb; mu++) {
            int A = orbs[mu].atom_idx;

            // Diagonal
            double f_diag = -orbs[mu].half_IP;
            f_diag += (P_AA[A] - params[A].z_valence
                       - (P_alpha(mu,mu) - 0.5)) * gamma(A, A);
            for (int B = 0; B < n_atoms; B++) {
                if (B != A)
                    f_diag += (P_AA[B] - params[B].z_valence) * gamma(A, B);
            }
            F(mu, mu) = f_diag;

            // Off-diagonal
            for (int nu = mu+1; nu < n_orb; nu++) {
                int B = orbs[nu].atom_idx;
                double f_off;
                if (A == B) {
                    // Same atom: ZDO → 0
                    f_off = 0.0;
                } else {
                    f_off = -0.5*(beta_orb[mu]+beta_orb[nu])*S(mu,nu)
                            - P_alpha(mu,nu)*gamma(A,B);
                }
                F(mu,nu) = f_off; F(nu,mu) = f_off;
            }
        }

        // Diagonalize F (ZDO: regular eigenvalue problem, S=I)
        // F * C = C * epsilon  (Lec 13, slide 8-9)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("cndo_energy: Fock diagonalization failed.");

        Eigen::MatrixXd C = solver.eigenvectors();

        // Update density matrix: P^total = 2 * C_occ * C_occ^T
        // P^alpha = C_occ * C_occ^T  (Lec 13, slide 8)
        Eigen::MatrixXd C_occ = C.leftCols(n_occ);
        Eigen::MatrixXd P_tot_new = 2.0 * C_occ * C_occ.transpose();

        // Total electronic energy: E = ½ Tr[P^tot * (H^core + F)]
        double E_elec = 0.5 * (P_tot_new.cwiseProduct(H_core + F)).sum();
        double E_total = E_elec + E_nuc;

        // Convergence check on both energy and density matrix
        double delta_E   = std::abs(E_total - E_old);
        double delta_dens = (P_tot_new - P_tot).norm() / n_orb;

        P_tot = P_tot_new;
        E_old = E_total;

        if (delta_E < CONV_ENER && delta_dens < CONV_DENS) {
            return E_total;
        }
    }

    std::cerr << "[WARNING] cndo_energy: SCF not converged after "
              << MAX_ITER << " iterations. Returning best estimate.\n";
    return E_old;
}

// ─── main ────────────────────────────────────────────────────────────────────
// Reads one .xyz file, runs CNDO/2 SCF, prints total energy (Hartree) to stdout.
int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: ./cndo_energy molecule.xyz\n";
            return 1;
        }
        Molecule mol;
        mol.read_xyz(argv[1]);
        double E = cndo2_total_energy(mol);
        std::cout << E << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
