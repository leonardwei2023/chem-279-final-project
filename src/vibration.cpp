#include "vibration.h"
#include "molecule.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static double get_atomic_mass(const std::string& atom) {
    if (atom == "H") return 1.00784;
    if (atom == "C") return 12.011;
    if (atom == "N") return 14.007;
    if (atom == "O") return 15.999;
    if (atom == "F") return 18.998;
    if (atom == "Cl" || atom == "CL") return 35.45;

    std::cerr << "Warning: unknown atom " << atom
              << ", using mass = 1.0\n";
    return 1.0;
}

static Eigen::MatrixXd read_hessian(const std::string& filename, int size) {
    Eigen::MatrixXd H(size, size);
    std::ifstream fin(filename);

    if (!fin) {
        throw std::runtime_error("Could not open Hessian file: " + filename);
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!(fin >> H(i, j))) {
                throw std::runtime_error("Invalid Hessian file format.");
            }
        }
    }

    return H;
}

static Eigen::MatrixXd mass_weight_hessian(
    const Eigen::MatrixXd& H,
    const Molecule& mol
) {
    int n_atoms = mol.atoms.size();
    int dim = 3 * n_atoms;

    Eigen::MatrixXd MW(dim, dim);

    for (int i = 0; i < dim; i++) {
        int ai = i / 3;
        double mi = get_atomic_mass(mol.atoms[ai].symbol);

        for (int j = 0; j < dim; j++) {
            int aj = j / 3;
            double mj = get_atomic_mass(mol.atoms[aj].symbol);

            MW(i, j) = H(i, j) / std::sqrt(mi * mj);
        }
    }

    return MW;
}

std::vector<double> compute_vibrational_frequencies(
    const Molecule& mol,
    const std::string& hessian_file
) {
    int n_atoms = mol.atoms.size();
    int dim = 3 * n_atoms;

    Eigen::MatrixXd H = read_hessian(hessian_file, dim);

    // Symmetrize
    H = 0.5 * (H + H.transpose());

    // Mass weighting
    Eigen::MatrixXd MW = mass_weight_hessian(H, mol);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(MW);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed.");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    std::vector<double> all_freqs;

    // 🔥 FINAL CORRECT CONVERSION
    const double AU_TO_CM = 5140.48;

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);

        // remove translation/rotation noise
        if (std::abs(lambda) < 1e-10) continue;

        double freq;

        if (lambda > 0.0) {
            freq = std::sqrt(lambda) * AU_TO_CM;
        } else {
            freq = -std::sqrt(std::abs(lambda)) * AU_TO_CM;
        }

        all_freqs.push_back(freq);
    }

    std::sort(all_freqs.begin(), all_freqs.end());

    int expected_modes;
    if (n_atoms == 2) {
        expected_modes = 1;
    } else {
        expected_modes = 3 * n_atoms - 6;
    }

    std::vector<double> vibrational_modes;

    for (double f : all_freqs) {
        if (std::abs(f) > 1.0) {
            vibrational_modes.push_back(f);
        }
    }

    if ((int)vibrational_modes.size() > expected_modes) {
        vibrational_modes.erase(
            vibrational_modes.begin(),
            vibrational_modes.end() - expected_modes
        );
    }

    return vibrational_modes;
}

void print_vibrational_frequencies(
    const Molecule& mol,
    const std::string& xyz_file,
    const std::string& hessian_file
) {
    std::vector<double> freqs =
        compute_vibrational_frequencies(mol, hessian_file);

    int n_atoms = mol.atoms.size();
    int expected_modes = (n_atoms == 2) ? 1 : 3 * n_atoms - 6;

    std::cout << "==============================\n";
    std::cout << "Molecule: " << mol.name << "\n";
    std::cout << "File: " << xyz_file << "\n";
    std::cout << "Atoms: " << n_atoms << "\n";
    std::cout << "==============================\n";
    std::cout << "Expected vibrational modes: "
              << expected_modes << "\n";
    std::cout << "Vibrational frequencies (cm^-1):\n";

    for (size_t i = 0; i < freqs.size(); i++) {
        std::cout << "Mode " << i + 1 << ": "
                  << freqs[i] << " cm^-1\n";
    }
}
