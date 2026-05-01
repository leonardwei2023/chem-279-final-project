#include "vibration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

bool Vibrations::is_linear_molecule(const Molecule& molecule) const {
    return molecule.get_num_atoms() == 2;
}

void Vibrations::compute(const Molecule& molecule, const Hessian& hessian) {
    frequencies.clear();

    int n_atoms = molecule.get_num_atoms();
    int dim = molecule.get_num_coordinates();

    if (hessian.get_size() != dim) {
        throw std::runtime_error("Hessian size does not match molecule coordinates.");
    }

    std::vector<double> masses = molecule.get_masses();

    Eigen::MatrixXd H(dim, dim);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            H(i, j) = hessian.get_value(i, j);
        }
    }

    H = 0.5 * (H + H.transpose());

    Eigen::MatrixXd MW(dim, dim);

    for (int i = 0; i < dim; i++) {
        int atom_i = i / 3;
        double mi = masses[atom_i];

        for (int j = 0; j < dim; j++) {
            int atom_j = j / 3;
            double mj = masses[atom_j];

            MW(i, j) = H(i, j) / std::sqrt(mi * mj);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(MW);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed.");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    std::vector<double> all_freqs;

    const double AU_TO_CM = 5140.48;

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);

        if (std::abs(lambda) < 1.0e-10) {
            continue;
        }

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

    if (is_linear_molecule(molecule)) {
        expected_modes = 3 * n_atoms - 5;
    } else {
        expected_modes = 3 * n_atoms - 6;
    }

    for (double f : all_freqs) {
        if (std::abs(f) > 1.0) {
            frequencies.push_back(f);
        }
    }

    if (static_cast<int>(frequencies.size()) > expected_modes) {
        frequencies.erase(
            frequencies.begin(),
            frequencies.end() - expected_modes
        );
    }
}

void Vibrations::print_frequencies() const {
    std::cout << "Vibrational frequencies (cm^-1):\n";

    for (size_t i = 0; i < frequencies.size(); i++) {
        std::cout << "Mode " << i + 1 << ": "
                  << frequencies[i] << " cm^-1\n";
    }
}

void Vibrations::write_frequencies(const std::string& filename) const {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not write frequency file: " + filename);
    }

    for (double f : frequencies) {
        file << f << "\n";
    }
}

std::vector<double> Vibrations::get_frequencies() const {
    return frequencies;
}
