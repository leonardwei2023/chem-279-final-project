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
    int n_atoms = static_cast<int>(mol.atoms.size());
    int dim = 3 * n_atoms;

    Eigen::MatrixXd MW(dim, dim);

    for (int i = 0; i < dim; i++) {
        int atom_i = i / 3;
        double mi = get_atomic_mass(mol.atoms[atom_i].symbol);

        for (int j = 0; j < dim; j++) {
            int atom_j = j / 3;
            double mj = get_atomic_mass(mol.atoms[atom_j].symbol);

            MW(i, j) = H(i, j) / std::sqrt(mi * mj);
        }
    }

    return MW;
}

std::vector<double> compute_vibrational_frequencies(
    const Molecule& mol,
    const std::string& hessian_file
) {
    int n_atoms = static_cast<int>(mol.atoms.size());
    int dim = 3 * n_atoms;

    Eigen::MatrixXd H = read_hessian(hessian_file, dim);

    // Remove small numerical asymmetry.
    H = 0.5 * (H + H.transpose());

    // Mass-weight Hessian: H_ij / sqrt(m_i m_j)
    Eigen::MatrixXd MW = mass_weight_hessian(H, mol);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(MW);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed.");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    std::vector<double> all_freqs;

    /*
      Unit conversion:

      If Hessian is in Hartree / Bohr^2 after mass weighting by amu,
      then eigenvalue units are Eh / (amu * Bohr^2).

      omega^2 = lambda * Eh_to_J / (amu_to_kg * bohr_to_m^2)
      frequency(cm^-1) = omega / (2*pi*c)
    */
    const double Eh_to_J = 4.3597447222071e-18;
    const double bohr_to_m = 5.29177210903e-11;
    const double amu_to_kg = 1.66053906660e-27;
    const double c_cm_s = 2.99792458e10;
    const double two_pi = 2.0 * M_PI;

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);

        // Skip tiny translational/rotational numerical noise.
        if (std::abs(lambda) < 1.0e-12) {
            continue;
        }

        double omega2 =
            lambda * Eh_to_J / (amu_to_kg * bohr_to_m * bohr_to_m);

        double freq_cm = 0.0;

        if (omega2 >= 0.0) {
            freq_cm = std::sqrt(omega2) / (two_pi * c_cm_s);
        } else {
            freq_cm = -std::sqrt(std::abs(omega2)) / (two_pi * c_cm_s);
        }

        all_freqs.push_back(freq_cm);
    }

    std::sort(all_freqs.begin(), all_freqs.end());

    int expected_modes = 0;
    if (n_atoms == 1) {
        expected_modes = 0;
    } else if (n_atoms == 2) {
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

    // Keep the largest expected_modes frequencies.
    if (static_cast<int>(vibrational_modes.size()) > expected_modes) {
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

    int n_atoms = static_cast<int>(mol.atoms.size());

    int expected_modes = 0;
    if (n_atoms == 1) {
        expected_modes = 0;
    } else if (n_atoms == 2) {
        expected_modes = 1;
    } else {
        expected_modes = 3 * n_atoms - 6;
    }

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
