#include "vibration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

bool Vibrations::is_linear_molecule(const Molecule& molecule) const {
    int n = molecule.get_num_atoms();

    if (n <= 2) {
        return true;
    }

    std::vector<double> c = molecule.get_coordinates();

    double ax = c[3] - c[0];
    double ay = c[4] - c[1];
    double az = c[5] - c[2];

    double norm_a = std::sqrt(ax * ax + ay * ay + az * az);

    if (norm_a < 1.0e-12) {
        return false;
    }

    for (int atom = 2; atom < n; atom++) {
        double bx = c[3 * atom] - c[0];
        double by = c[3 * atom + 1] - c[1];
        double bz = c[3 * atom + 2] - c[2];

        double cx = ay * bz - az * by;
        double cy = az * bx - ax * bz;
        double cz = ax * by - ay * bx;

        double cross_norm = std::sqrt(cx * cx + cy * cy + cz * cz);

        if (cross_norm > 1.0e-6) {
            return false;
        }
    }

    return true;
}

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
        std::cout << "Eigenvalue calculation failed.\n";
        return;
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    std::vector<double> vals;

    for (int i = 0; i < eigenvalues.size(); i++) {
        vals.push_back(eigenvalues(i));
    }

    std::sort(vals.begin(), vals.end());

    int num_modes;

    if (is_linear_molecule(molecule)) {
        num_modes = 3 * num_atoms - 5;
    }
    else {
        num_modes = 3 * num_atoms - 6;
    }

    int start_index = vals.size() - num_modes;

    
    /*
       Unit Conversions:
       Hessian units = Hartree / Angstrom^2
       Mass units    = amu
       Output        = cm^-1
    */

    
    const double conversion_to_cm = 2721.14;

    for (int i = start_index; i < vals.size(); i++) {
        double lambda = vals[i];

        if (lambda > 1.0e-10) {
            double freq = std::sqrt(lambda) * conversion_to_cm;
            frequencies.push_back(freq);
        }
        else {
            frequencies.push_back(0.0);
        }
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

    for (double freq : frequencies) {
        file << freq << "\n";
    }
}

std::vector<double> Vibrations::get_frequencies() const {
    return frequencies;
}
