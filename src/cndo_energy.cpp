#include "molecule.h"

#include <Eigen/Dense>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct AO {
    int atom;
    std::string type;
};

struct Params {
    int valence;
    double beta_ev;
    double I_ev;
    double A_ev;
};

static const double EV_TO_HARTREE = 1.0 / 27.211386245988;
static const double ANG_TO_BOHR = 1.889726124565062;

Params get_params(const std::string& s) {
    if (s == "H")  return {1, -9.0, 13.60, 0.754};
    if (s == "C")  return {4, -21.0, 21.40, 1.260};
    if (s == "N")  return {5, -25.0, 26.00, 1.800};
    if (s == "O")  return {6, -31.0, 32.30, 1.460};
    if (s == "F")  return {7, -39.0, 40.00, 3.400};
    if (s == "Cl") return {7, -15.0, 25.23, 3.610};

    return {1, -9.0, 13.60, 0.754};
}

double distance_bohr(const std::vector<double>& c, int a, int b) {
    double dx = c[3 * a]     - c[3 * b];
    double dy = c[3 * a + 1] - c[3 * b + 1];
    double dz = c[3 * a + 2] - c[3 * b + 2];

    return std::sqrt(dx * dx + dy * dy + dz * dz) * ANG_TO_BOHR;
}

std::vector<AO> build_basis(const std::vector<std::string>& symbols) {
    std::vector<AO> basis;

    for (int a = 0; a < static_cast<int>(symbols.size()); a++) {
        if (symbols[a] == "H") {
            basis.push_back({a, "s"});
        } else {
            basis.push_back({a, "s"});
            basis.push_back({a, "px"});
            basis.push_back({a, "py"});
            basis.push_back({a, "pz"});
        }
    }

    return basis;
}

double overlap_approx(
    const AO& mu,
    const AO& nu,
    const std::vector<double>& coords
) {
    if (mu.atom == nu.atom) {
        return (mu.type == nu.type) ? 1.0 : 0.0;
    }

    double R = distance_bohr(coords, mu.atom, nu.atom);
    if (R < 1.0e-12) return 0.0;

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

double gamma_AB(
    const std::string& A,
    const std::string& B,
    double R_bohr
) {
    Params pA = get_params(A);
    Params pB = get_params(B);

    double etaA = 0.5 * (pA.I_ev - pA.A_ev) * EV_TO_HARTREE;
    double etaB = 0.5 * (pB.I_ev - pB.A_ev) * EV_TO_HARTREE;

    double gammaAA = std::max(etaA, 1.0e-8);
    double gammaBB = std::max(etaB, 1.0e-8);

    if (R_bohr < 1.0e-10) {
        return gammaAA;
    }

    double avg = 0.5 * (gammaAA + gammaBB);

    return 1.0 / std::sqrt(R_bohr * R_bohr + 1.0 / (avg * avg));
}

double nuclear_repulsion(
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    double e = 0.0;

    for (int A = 0; A < static_cast<int>(symbols.size()); A++) {
        double ZA = get_params(symbols[A]).valence;

        for (int B = A + 1; B < static_cast<int>(symbols.size()); B++) {
            double ZB = get_params(symbols[B]).valence;
            double R = distance_bohr(coords, A, B);

            if (R > 1.0e-12) {
                e += ZA * ZB / R;
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
    int nbf = static_cast<int>(basis.size());
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nbf, nbf);

    for (int mu = 0; mu < nbf; mu++) {
        int A = basis[mu].atom;
        Params pA = get_params(symbols[A]);

        H(mu, mu) = -0.5 * (pA.I_ev + pA.A_ev) * EV_TO_HARTREE;
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = mu + 1; nu < nbf; nu++) {
            int A = basis[mu].atom;
            int B = basis[nu].atom;

            if (A == B) continue;

            Params pA = get_params(symbols[A]);
            Params pB = get_params(symbols[B]);

            double betaA = pA.beta_ev * EV_TO_HARTREE;
            double betaB = pB.beta_ev * EV_TO_HARTREE;
            double Smunu = overlap_approx(basis[mu], basis[nu], coords);

            H(mu, nu) = -0.5 * (betaA + betaB) * Smunu;
            H(nu, mu) = H(mu, nu);
        }
    }

    return H;
}

Eigen::MatrixXd make_density(
    const Eigen::MatrixXd& C,
    int electrons
) {
    int nbf = C.rows();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nbf, nbf);

    int occupied = electrons / 2;

    for (int i = 0; i < occupied; i++) {
        P += 2.0 * C.col(i) * C.col(i).transpose();
    }

    return P;
}

Eigen::MatrixXd build_fock(
    const Eigen::MatrixXd& H,
    const Eigen::MatrixXd& Ptot,
    const std::vector<AO>& basis,
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords
) {
    int nbf = static_cast<int>(basis.size());
    int nat = static_cast<int>(symbols.size());

    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(nbf, nbf);

    std::vector<double> PAA(nat, 0.0);

    for (int mu = 0; mu < nbf; mu++) {
        int A = basis[mu].atom;
        PAA[A] += Ptot(mu, mu);
    }

    for (int mu = 0; mu < nbf; mu++) {
        int A = basis[mu].atom;
        Params pA = get_params(symbols[A]);

        double gammaAA = gamma_AB(symbols[A], symbols[A], 0.0);
        double P_alpha_mumu = 0.5 * Ptot(mu, mu);

        double value = -0.5 * (pA.I_ev + pA.A_ev) * EV_TO_HARTREE;

        value += ((PAA[A] - pA.valence) - (P_alpha_mumu - 0.5)) * gammaAA;

        for (int B = 0; B < nat; B++) {
            if (B == A) continue;

            double gammaAB = gamma_AB(
                symbols[A],
                symbols[B],
                distance_bohr(coords, A, B)
            );

            Params pB = get_params(symbols[B]);
            value += (PAA[B] - pB.valence) * gammaAB;
        }

        F(mu, mu) = value;
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = mu + 1; nu < nbf; nu++) {
            int A = basis[mu].atom;
            int B = basis[nu].atom;

            double value = 0.0;

            if (A != B) {
                Params pA = get_params(symbols[A]);
                Params pB = get_params(symbols[B]);

                double betaA = pA.beta_ev * EV_TO_HARTREE;
                double betaB = pB.beta_ev * EV_TO_HARTREE;
                double Smunu = overlap_approx(basis[mu], basis[nu], coords);

                double gammaAB = gamma_AB(
                    symbols[A],
                    symbols[B],
                    distance_bohr(coords, A, B)
                );

                double P_alpha_munu = 0.5 * Ptot(mu, nu);

                value = -0.5 * (betaA + betaB) * Smunu
                        - P_alpha_munu * gammaAB;
            }

            F(mu, nu) = value;
            F(nu, mu) = value;
        }
    }

    return F;
}

double electronic_energy(
    const Eigen::MatrixXd& P,
    const Eigen::MatrixXd& H,
    const Eigen::MatrixXd& F
) {
    double e = 0.0;

    for (int mu = 0; mu < P.rows(); mu++) {
        for (int nu = 0; nu < P.cols(); nu++) {
            e += 0.5 * P(mu, nu) * (H(mu, nu) + F(mu, nu));
        }
    }

    return e;
}

double compute_cndo2_energy(const Molecule& mol) {
    std::vector<std::string> symbols = mol.get_symbols();
    std::vector<double> coords = mol.get_coordinates();

    std::vector<AO> basis = build_basis(symbols);

    int nbf = static_cast<int>(basis.size());
    int electrons = 0;

    for (const auto& s : symbols) {
        electrons += get_params(s).valence;
    }

    if (electrons % 2 != 0) {
        throw std::runtime_error(
            "This implementation supports closed-shell molecules only."
        );
    }

    Eigen::MatrixXd H = build_core_hamiltonian(basis, symbols, coords);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> guess_solver(H);
    if (guess_solver.info() != Eigen::Success) {
        throw std::runtime_error("Initial Huckel guess diagonalization failed.");
    }

    Eigen::MatrixXd P = make_density(guess_solver.eigenvectors(), electrons);

    double old_energy = 1.0e100;
    double total_energy = 0.0;

    const int max_iter = 200;
    const double tolerance = 1.0e-9;
    const double damping = 0.50;

    for (int iter = 0; iter < max_iter; iter++) {
        Eigen::MatrixXd F = build_fock(H, P, basis, symbols, coords);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("SCF Fock diagonalization failed.");
        }

        Eigen::MatrixXd P_new = make_density(solver.eigenvectors(), electrons);

        P = damping * P_new + (1.0 - damping) * P;

        double e_elec = electronic_energy(P, H, F);
        double e_nuc = nuclear_repulsion(symbols, coords);

        total_energy = e_elec + e_nuc;

        if (std::abs(total_energy - old_energy) < tolerance) {
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

        Molecule mol;
        mol.read_xyz(argv[1]);

        double energy = compute_cndo2_energy(mol);

        std::cout << energy << std::endl;
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
