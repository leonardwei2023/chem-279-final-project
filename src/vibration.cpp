
#include "vibrations.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

void Vibrations::compute(const Molecule& molecule, const Hessian& hessian) {
    int num_atoms = molecule.get_num_atoms();
    int size = 3 * num_atoms;

    std::vector<double> masses = molecule.get_masses();
    std::vector<std::vector<double>> H = hessian.get_matrix();

    Eigen::MatrixXd mass_weighted(size, size);

    for (int i = 0; i < size; i++) {
        int atom_i = i / 3;

        for (int j = 0; j < size; j++) {
            int atom_j = j / 3;

            double mass_factor = std::sqrt(masses[atom_i] * masses[atom_j]);
            mass_weighted(i, j) = H[i][j] / mass_factor;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mass_weighted);

    frequencies.clear();

    if (solver.info() != Eigen::Success) {
        std::cout << "Eigenvalue calculation failed." << std::endl;
        return;
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);

        if (lambda > 1.0e-8) {
            double freq = std::sqrt(lambda);
            frequencies.push_back(freq);
        }
    }
}

void Vibrations::print_frequencies() const {
    std::cout << "Vibrational frequencies:" << std::endl;

    for (size_t i = 0; i < frequencies.size(); i++) {
        std::cout << "Mode " << i + 1 << ": " << frequencies[i] << std::endl;
    }
}
