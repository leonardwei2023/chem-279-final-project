#include "molecule.h"

#include <Eigen/Dense>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct AO {
    int atom;
    std::string type;
};

struct Params {
    int valence;
    double zeta;
    double beta_eV;
    double ip_s_eV;
    double ip_p_eV;
    double ea_s_eV;
    double ea_p_eV;
};

static const double EV_TO_HARTREE = 1.0 / 27.211386245988;
static const double ANG_TO_BOHR = 1.889726124565062;

Params get_params(const std::string& s) {
    if (s == "H")  return {1, 1.20, -9.0, 13.60, 13.60, 0.754, 0.754};
    if (s == "C")  return {4, 1.625, -21.0, 21.40, 11.40, 1.260, 1.260};
    if (s == "N")  return {5, 1.950, -25.0, 26.00, 13.40, 1.800, 1.800};
    if (s == "O")  return {6, 2.275, -31.0, 32.30, 14.80, 1.460, 1.460};
    if (s == "F")  return {7, 2.550, -39.0, 40.00, 18.10, 3.400, 3.400};
    if (s == "Cl") return {7, 2.033, -15.0, 25.23, 15.03, 3.610, 3.610};

    return {1, 1.20, -9.0, 13.60, 13.60, 0.754, 0.754};
}

double distance_bohr(const std::vector<double>& c, int a, int b) {
    double dx = c[3 * a]     - c[3 * b];
    double dy = c[3 * a + 1] - c[3 * b + 1];
    double dz = c[3 * a + 2] - c[3 * b + 2];

    return std::sqrt(dx * dx + dy * dy + dz * dz) * ANG_TO_BOHR;
}

std::vector<AO> build_basis(const std::vector<std::string>& symbols) {
    std::vector<AO> basis;

    for (int i = 0; i < (int)symbols.size(); i++) {
        if (symbols[i] == "H") {
            basis.push_back({i, "s"});
        } else {
            basis.push_back({i, "s"});
            basis.push_back({i, "px"});
            basis.push_back({i, "py"});
            basis.push_back({i, "pz"});
        }
    }

    return basis;
}

double gamma_AB(const std::string& A, const std::string& B, double R_bohr) {
    Params pA = get_params(A);
    Params pB = get_params(B);

    double etaA = 0.5 * ((pA.ip_s_eV - pA.ea_s_eV) * EV_TO_HARTREE);
    double etaB = 0.5 * ((pB.ip_s_eV - pB.ea_s_eV) * EV_TO_HARTREE);

    double aA = 1.0 / std::max(etaA, 1.0e-6);
    double aB = 1.0 / std::max(etaB, 1.0e-6);

    return 1.0 / std::sqrt(R_bohr * R_bohr + 0.25 * (aA + aB) * (aA + aB));
}

double overlap_approx(
    const AO& mu,
    const AO& nu,
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    if (mu.atom == nu.atom) {
        return (mu.type == nu.type) ? 1.0 : 0.0;
    }

    double R = distance_bohr(coords, mu.atom, nu.atom);
    if (R < 1.0e-8) return 0.0;

    double S = std::exp(-0.65 * R);

    double dx = (coords[3 * nu.atom]     - coords[3 * mu.atom]) * ANG_TO_BOHR / R;
    double dy = (coords[3 * nu.atom + 1] - coords[3 * mu.atom + 1]) * ANG_TO_BOHR / R;
    double dz = (coords[3 * nu.atom + 2] - coords[3 * mu.atom + 2]) * ANG_TO_BOHR / R;

    auto orient = [&](const std::string& t) {
        if (t == "px") return dx;
        if (t == "py") return dy;
        if (t == "pz") return dz;
        return 1.0;
    };

    return S * orient(mu.type) * orient(nu.type);
}

double core_diagonal(const AO& ao, const std::vector<std::string>& symbols) {
    Params p = get_params(symbols[ao.atom]);

    if (ao.type == "s") {
        return -p.ip_s_eV * EV_TO_HARTREE;
    }

    return -p.ip_p_eV * EV_TO_HARTREE;
}

double nuclear_repulsion(
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    double e = 0.0;

    for (int i = 0; i < (int)symbols.size(); i++) {
        int Zi = get_params(symbols[i]).valence;

        for (int j = i + 1; j < (int)symbols.size(); j++) {
            int Zj = get_params(symbols[j]).valence;
            double R = distance_bohr(coords, i, j);

            if (R > 1.0e-8) {
                e += Zi * Zj / R;
            }
        }
    }

    return e;
}

Eigen::MatrixXd build_core_hamiltonian(
    const std::vector<AO>& basis,
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    int nbf = basis.size();
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nbf, nbf);

    for (int mu = 0; mu < nbf; mu++) {
        H(mu, mu) = core_diagonal(basis[mu], symbols);
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = mu + 1; nu < nbf; nu++) {
            int A = basis[mu].atom;
            int B = basis[nu].atom;

            if (A == B) continue;

            Params pA = get_params(symbols[A]);
            Params pB = get_params(symbols[B]);

            double beta = 0.5 * (pA.beta_eV + pB.beta_eV) * EV_TO_HARTREE;
            double S = overlap_approx(basis[mu], basis[nu], symbols, coords);

            H(mu, nu) = beta * S;
            H(nu, mu) = H(mu, nu);
        }
    }

    return H;
}

Eigen::MatrixXd build_fock(
    const Eigen::MatrixXd& H,
    const Eigen::MatrixXd& P,
    const std::vector<AO>& basis,
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    int nbf = basis.size();
    int nat = symbols.size();

    std::vector<double> P_atom(nat, 0.0);

    for (int mu = 0; mu < nbf; mu++) {
        P_atom[basis[mu].atom] += P(mu, mu);
    }

    Eigen::MatrixXd F = H;

    for (int mu = 0; mu < nbf; mu++) {
        int A = basis[mu].atom;

        double val = H(mu, mu);

        for (int B = 0; B < nat; B++) {
            double R = (A == B) ? 0.0 : distance_bohr(coords, A, B);
            double gamma = gamma_AB(symbols[A], symbols[B], R);
            double ZB = get_params(symbols[B]).valence;

            val += (P_atom[B] - ZB) * gamma;
        }

        double gammaAA = gamma_AB(symbols[A], symbols[A], 0.0);
        val -= 0.5 * P(mu, mu) * gammaAA;

        F(mu, mu) = val;
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = mu + 1; nu < nbf; nu++) {
            int A = basis[mu].atom;
            int B = basis[nu].atom;

            double R = (A == B) ? 0.0 : distance_bohr(coords, A, B);
            double gamma = gamma_AB(symbols[A], symbols[B], R);

            F(mu, nu) = H(mu, nu) - 0.5 * P(mu, nu) * gamma;
            F(nu, mu) = F(mu, nu);
        }
    }

    return F;
}

Eigen::MatrixXd density_from_coefficients(
    const Eigen::MatrixXd& C,
    int electrons
) {
    int nbf = C.rows();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nbf, nbf);

    int doubly_occ = electrons / 2;
    bool odd = (electrons % 2 == 1);

    for (int i = 0; i < doubly_occ; i++) {
        P += 2.0 * C.col(i) * C.col(i).transpose();
    }

    if (odd && doubly_occ < C.cols()) {
        P += C.col(doubly_occ) * C.col(doubly_occ).transpose();
    }

    return P;
}

double electronic_energy(
    const Eigen::MatrixXd& P,
    const Eigen::MatrixXd& H,
    const Eigen::MatrixXd& F
) {
    double E = 0.0;

    for (int mu = 0; mu < P.rows(); mu++) {
        for (int nu = 0; nu < P.cols(); nu++) {
            E += 0.5 * P(mu, nu) * (H(mu, nu) + F(mu, nu));
        }
    }

    return E;
}

double compute_cndo2_energy(const Molecule& molecule) {
    std::vector<std::string> symbols = molecule.get_symbols();
    std::vector<double> coords = molecule.get_coordinates();

    std::vector<AO> basis = build_basis(symbols);

    int nbf = basis.size();
    int electrons = 0;

    for (const auto& s : symbols) {
        electrons += get_params(s).valence;
    }

    Eigen::MatrixXd H = build_core_hamiltonian(basis, symbols, coords);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nbf, nbf);
    Eigen::MatrixXd F = H;

    double old_energy = 1.0e100;
    double total_energy = 0.0;

    const int max_iter = 200;
    const double tol = 1.0e-9;
    const double damping = 0.35;

    for (int iter = 0; iter < max_iter; iter++) {
        F = build_fock(H, P, basis, symbols, coords);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("CNDO/2 SCF diagonalization failed.");
        }

        Eigen::MatrixXd C = solver.eigenvectors();
        Eigen::MatrixXd P_new = density_from_coefficients(C, electrons);

        P = damping * P_new + (1.0 - damping) * P;

        double E_elec = electronic_energy(P, H, F);
        double E_nuc = nuclear_repulsion(symbols, coords);
        total_energy = E_elec + E_nuc;

        if (std::abs(total_energy - old_energy) < tol) {
            break;
        }

        old_energy = total_energy;
    }

    return total_energy;
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cout << "Usage: ./cndo_energy molecule.xyz\n";
            return 1;
        }

        Molecule molecule;
        molecule.read_xyz(argv[1]);

        double energy = compute_cndo2_energy(molecule);

        std::cout << energy << std::endl;
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
