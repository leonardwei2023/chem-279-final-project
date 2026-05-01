#include "molecule.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static constexpr double EV_TO_HA    = 0.03674930;
static constexpr double ANG_TO_BOHR = 1.8897261245;
static constexpr double PI          = 3.14159265358979323846;

struct GTO {
    double alpha;
    double coeff;
};

static const GTO STO3G_1s[3] = {
    {2.2277, 0.1543}, {0.4058, 0.5353}, {0.1098, 0.4446}
};

static const GTO STO3G_2s[3] = {
    {2.5820, -0.0999}, {0.1567, 0.3995}, {0.0611, 0.7002}
};

static const GTO STO3G_2p[3] = {
    {2.5820, 0.1559}, {0.1567, 0.6077}, {0.0611, 0.3920}
};

struct AtomParams {
    int z_valence;
    int n_orbitals;
    double I_s;
    double A_s;
    double I_p;
    double A_p;
    double beta;
    double zeta_s;
    double zeta_p;
};

static AtomParams get_params(const std::string& sym) {
    if (sym == "H")  return {1, 1, 13.06,  3.69,  0.0,   0.0,   -9.0,  1.200, 0.000};
    if (sym == "C")  return {4, 4, 14.051, 5.572, 5.572, 2.275, -21.0, 1.625, 1.625};
    if (sym == "N")  return {5, 4, 19.316, 7.275, 7.275, 2.862, -25.0, 1.950, 1.950};
    if (sym == "O")  return {6, 4, 25.390, 9.111, 9.111, 3.524, -31.0, 2.275, 2.275};
    if (sym == "F")  return {7, 4, 32.272, 11.08, 11.08, 4.270, -39.0, 2.425, 2.425};
    if (sym == "Cl") return {7, 4, 24.350, 3.890, 3.890, 1.430, -15.0, 3.500, 2.033};

    throw std::runtime_error("cndo_energy: unsupported element: " + sym);
}

enum OrbType {
    ORB_1S,
    ORB_2S,
    ORB_2PX,
    ORB_2PY,
    ORB_2PZ
};

struct OrbInfo {
    int atom_idx;
    OrbType type;
    double zeta;
    double half_IP;
};

static std::vector<OrbInfo> make_orbitals(
    const Molecule& mol,
    const std::vector<AtomParams>& params
) {
    std::vector<OrbInfo> orbs;
    int n = mol.get_num_atoms();

    for (int A = 0; A < n; A++) {
        const AtomParams& p = params[A];

        if (p.n_orbitals == 1) {
            double half_s = 0.5 * (p.I_s + p.A_s) * EV_TO_HA;
            orbs.push_back({A, ORB_1S, p.zeta_s, half_s});
        } else {
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

static double prim_ss(
    double a,
    double b,
    const double Ra[3],
    const double Rb[3]
) {
    double dx = Ra[0] - Rb[0];
    double dy = Ra[1] - Rb[1];
    double dz = Ra[2] - Rb[2];

    double R2 = dx * dx + dy * dy + dz * dz;

    return std::pow(PI / (a + b), 1.5) *
           std::exp(-a * b / (a + b) * R2);
}

static double prim_ps(
    int pi,
    double a,
    double b,
    const double Ra[3],
    const double Rb[3]
) {
    double ab = a + b;

    double P[3];
    for (int k = 0; k < 3; k++) {
        P[k] = (a * Ra[k] + b * Rb[k]) / ab;
    }

    double dR2 = 0.0;
    for (int k = 0; k < 3; k++) {
        double d = Ra[k] - Rb[k];
        dR2 += d * d;
    }

    double PmA = P[pi] - Ra[pi];

    return PmA * std::pow(PI / ab, 1.5) *
           std::exp(-a * b / ab * dR2);
}

static double prim_pp(
    int pi,
    int pj,
    double a,
    double b,
    const double Ra[3],
    const double Rb[3]
) {
    double ab = a + b;

    double P[3];
    for (int k = 0; k < 3; k++) {
        P[k] = (a * Ra[k] + b * Rb[k]) / ab;
    }

    double dR2 = 0.0;
    for (int k = 0; k < 3; k++) {
        double d = Ra[k] - Rb[k];
        dR2 += d * d;
    }

    double PmA = P[pi] - Ra[pi];
    double PmB = P[pj] - Rb[pj];

    double s = std::pow(PI / ab, 1.5) *
               std::exp(-a * b / ab * dR2);

    return (PmA * PmB + (pi == pj ? 0.5 / ab : 0.0)) * s;
}

static double sto3g_overlap(
    const OrbInfo& mu,
    const double* Rmu,
    const OrbInfo& nu,
    const double* Rnu
) {
    auto gtos = [](OrbType t) -> const GTO* {
        if (t == ORB_1S) return STO3G_1s;
        if (t == ORB_2S) return STO3G_2s;
        return STO3G_2p;
    };

    const GTO* ga = gtos(mu.type);
    const GTO* gb = gtos(nu.type);

    bool is_p_a = (mu.type == ORB_2PX ||
                   mu.type == ORB_2PY ||
                   mu.type == ORB_2PZ);

    bool is_p_b = (nu.type == ORB_2PX ||
                   nu.type == ORB_2PY ||
                   nu.type == ORB_2PZ);

    int pi_a = (mu.type == ORB_2PX) ? 0 :
               (mu.type == ORB_2PY) ? 1 : 2;

    int pi_b = (nu.type == ORB_2PX) ? 0 :
               (nu.type == ORB_2PY) ? 1 : 2;

    double S = 0.0;

    for (int j = 0; j < 3; j++) {
        double alpha_a = ga[j].alpha * mu.zeta * mu.zeta;

        for (int k = 0; k < 3; k++) {
            double alpha_b = gb[k].alpha * nu.zeta * nu.zeta;

            double s_prim = 0.0;

            if (!is_p_a && !is_p_b) {
                s_prim = prim_ss(alpha_a, alpha_b, Rmu, Rnu);
            } else if (is_p_a && !is_p_b) {
                s_prim = prim_ps(pi_a, alpha_a, alpha_b, Rmu, Rnu);
            } else if (!is_p_a && is_p_b) {
                s_prim = prim_ps(pi_b, alpha_b, alpha_a, Rnu, Rmu);
            } else {
                s_prim = prim_pp(pi_a, pi_b, alpha_a, alpha_b, Rmu, Rnu);
            }

            S += ga[j].coeff * gb[k].coeff * s_prim;
        }
    }

    return S;
}

static double gamma_onsite(const AtomParams& p) {
    return (p.I_s - p.A_s) * EV_TO_HA;
}

static double gamma_two_center(double R_bohr, double gAA, double gBB) {
    return 1.0 / (R_bohr + 2.0 / (gAA + gBB));
}

static void write_p_diagonal_file(
    const std::string& filename,
    int n_atoms,
    int n_orb,
    const std::vector<OrbInfo>& orbs,
    const Eigen::MatrixXd& P_tot
) {
    std::ofstream pfile(filename);

    if (!pfile.is_open()) {
        throw std::runtime_error("cndo_energy: could not write " + filename);
    }

    for (int A = 0; A < n_atoms; A++) {
        double population = 0.0;

        for (int mu = 0; mu < n_orb; mu++) {
            if (orbs[mu].atom_idx == A) {
                population += P_tot(mu, mu);
            }
        }

        pfile << population << "\n";
    }
}

double cndo2_total_energy(const Molecule& mol) {
    int n_atoms = mol.get_num_atoms();
    auto syms = mol.get_symbols();
    auto coords = mol.get_coordinates();

    std::vector<AtomParams> params(n_atoms);
    int n_elec = 0;

    for (int A = 0; A < n_atoms; A++) {
        params[A] = get_params(syms[A]);
        n_elec += params[A].z_valence;
    }

    if (n_elec % 2 != 0) {
        throw std::runtime_error(
            "cndo_energy: odd electron count — open shell not supported."
        );
    }

    int n_occ = n_elec / 2;

    std::vector<std::array<double, 3>> R(n_atoms);

    for (int A = 0; A < n_atoms; A++) {
        R[A][0] = coords[3 * A]     * ANG_TO_BOHR;
        R[A][1] = coords[3 * A + 1] * ANG_TO_BOHR;
        R[A][2] = coords[3 * A + 2] * ANG_TO_BOHR;
    }

    auto orbs = make_orbitals(mol, params);
    int n_orb = static_cast<int>(orbs.size());

    std::vector<double> beta_orb(n_orb);

    for (int mu = 0; mu < n_orb; mu++) {
        beta_orb[mu] = params[orbs[mu].atom_idx].beta * EV_TO_HA;
    }

    std::vector<double> gAA(n_atoms);

    for (int A = 0; A < n_atoms; A++) {
        gAA[A] = gamma_onsite(params[A]);
    }

    Eigen::MatrixXd gamma(n_atoms, n_atoms);

    for (int A = 0; A < n_atoms; A++) {
        gamma(A, A) = gAA[A];

        for (int B = A + 1; B < n_atoms; B++) {
            double dx = R[A][0] - R[B][0];
            double dy = R[A][1] - R[B][1];
            double dz = R[A][2] - R[B][2];

            double Rab = std::sqrt(dx * dx + dy * dy + dz * dz);
            double g = gamma_two_center(Rab, gAA[A], gAA[B]);

            gamma(A, B) = g;
            gamma(B, A) = g;
        }
    }

    Eigen::MatrixXd S = Eigen::MatrixXd::Identity(n_orb, n_orb);

    for (int mu = 0; mu < n_orb; mu++) {
        int A = orbs[mu].atom_idx;

        for (int nu = mu + 1; nu < n_orb; nu++) {
            int B = orbs[nu].atom_idx;

            double s = sto3g_overlap(
                orbs[mu],
                R[A].data(),
                orbs[nu],
                R[B].data()
            );

            S(mu, nu) = s;
            S(nu, mu) = s;
        }
    }

    Eigen::MatrixXd H_core = Eigen::MatrixXd::Zero(n_orb, n_orb);

    for (int mu = 0; mu < n_orb; mu++) {
        int A = orbs[mu].atom_idx;

        double h = -orbs[mu].half_IP;
        h += (params[A].z_valence - 0.5) * gamma(A, A);

        for (int B = 0; B < n_atoms; B++) {
            if (B != A) {
                h -= params[B].z_valence * gamma(A, B);
            }
        }

        H_core(mu, mu) = h;

        for (int nu = mu + 1; nu < n_orb; nu++) {
            int B = orbs[nu].atom_idx;

            if (A != B) {
                double hoff = -0.5 *
                              (beta_orb[mu] + beta_orb[nu]) *
                              S(mu, nu);

                H_core(mu, nu) = hoff;
                H_core(nu, mu) = hoff;
            }
        }
    }

    double E_nuc = 0.0;

    for (int A = 0; A < n_atoms; A++) {
        for (int B = A + 1; B < n_atoms; B++) {
            E_nuc += params[A].z_valence *
                     params[B].z_valence *
                     gamma(A, B);
        }
    }

    Eigen::MatrixXd P_tot = Eigen::MatrixXd::Zero(n_orb, n_orb);

    const int MAX_ITER = 300;
    const double CONV_DENS = 1.0e-8;
    const double CONV_ENER = 1.0e-10;

    double E_old = 1.0e30;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        std::vector<double> P_AA(n_atoms, 0.0);

        for (int mu = 0; mu < n_orb; mu++) {
            P_AA[orbs[mu].atom_idx] += P_tot(mu, mu);
        }

        Eigen::MatrixXd P_alpha = 0.5 * P_tot;
        Eigen::MatrixXd F(n_orb, n_orb);

        for (int mu = 0; mu < n_orb; mu++) {
            int A = orbs[mu].atom_idx;

            double f_diag = -orbs[mu].half_IP;

            f_diag += (
                P_AA[A] -
                params[A].z_valence -
                (P_alpha(mu, mu) - 0.5)
            ) * gamma(A, A);

            for (int B = 0; B < n_atoms; B++) {
                if (B != A) {
                    f_diag += (
                        P_AA[B] -
                        params[B].z_valence
                    ) * gamma(A, B);
                }
            }

            F(mu, mu) = f_diag;

            for (int nu = mu + 1; nu < n_orb; nu++) {
                int B = orbs[nu].atom_idx;

                double f_off = 0.0;

                if (A == B) {
                    f_off = 0.0;
                } else {
                    f_off = -0.5 *
                            (beta_orb[mu] + beta_orb[nu]) *
                            S(mu, nu)
                            - P_alpha(mu, nu) * gamma(A, B);
                }

                F(mu, nu) = f_off;
                F(nu, mu) = f_off;
            }
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("cndo_energy: Fock diagonalization failed.");
        }

        Eigen::MatrixXd C = solver.eigenvectors();
        Eigen::MatrixXd C_occ = C.leftCols(n_occ);
        Eigen::MatrixXd P_tot_new = 2.0 * C_occ * C_occ.transpose();

        double E_elec =
            0.5 * (P_tot_new.cwiseProduct(H_core + F)).sum();

        double E_total = E_elec + E_nuc;

        double delta_E = std::abs(E_total - E_old);
        double delta_dens = (P_tot_new - P_tot).norm() / n_orb;

        P_tot = P_tot_new;
        E_old = E_total;

        if (delta_E < CONV_ENER && delta_dens < CONV_DENS) {
            write_p_diagonal_file(
                "p_diagonal.dat",
                n_atoms,
                n_orb,
                orbs,
                P_tot
            );

            return E_total;
        }
    }

    std::cerr << "[WARNING] cndo_energy: SCF not converged after "
              << MAX_ITER << " iterations. Returning best estimate.\n";

    write_p_diagonal_file(
        "p_diagonal.dat",
        n_atoms,
        n_orb,
        orbs,
        P_tot
    );

    return E_old;
}

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
