#include "dipole.h"
#include "molecule.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int DipoleMoment::get_valence_electrons(const std::string& symbol) {
    if (symbol == "H")  return 1;
    if (symbol == "C")  return 4;
    if (symbol == "N")  return 5;
    if (symbol == "O")  return 6;
    if (symbol == "F")  return 7;
    if (symbol == "Cl") return 7;

    throw std::runtime_error("DipoleMoment: unknown element symbol: " + symbol);
}

void DipoleMoment::compute(
    const Molecule& molecule,
    const std::vector<double>& p_diagonal
) {
    int n_atoms = molecule.get_num_atoms();

    std::vector<std::string> symbols = molecule.get_symbols();
    std::vector<double> coords = molecule.get_coordinates();

    if (p_diagonal.empty()) {
        throw std::runtime_error(
            "DipoleMoment::compute: no SCF density provided. "
            "Run CNDO/2 SCF first to generate p_diagonal.dat."
        );
    }

    if (static_cast<int>(p_diagonal.size()) != n_atoms) {
        throw std::runtime_error(
            "DipoleMoment::compute: p_diagonal size does not match atom count."
        );
    }

    mu_debye = {0.0, 0.0, 0.0};

    for (int A = 0; A < n_atoms; A++) {
        double Z_A = static_cast<double>(get_valence_electrons(symbols[A]));
        double electron_population = p_diagonal[A];

        double net_charge = Z_A - electron_population;

        double x = coords[3 * A];
        double y = coords[3 * A + 1];
        double z = coords[3 * A + 2];

        mu_debye[0] += net_charge * x;
        mu_debye[1] += net_charge * y;
        mu_debye[2] += net_charge * z;
    }

    mu_debye[0] *= EA_TO_DEBYE;
    mu_debye[1] *= EA_TO_DEBYE;
    mu_debye[2] *= EA_TO_DEBYE;

    magnitude_debye = std::sqrt(
        mu_debye[0] * mu_debye[0] +
        mu_debye[1] * mu_debye[1] +
        mu_debye[2] * mu_debye[2]
    );
}

void DipoleMoment::print() const {
    std::cout << "\n--- Dipole Moment (CNDO/2 SCF Density) ---\n";
    std::cout << "  mu_x = " << mu_debye[0] << " Debye\n";
    std::cout << "  mu_y = " << mu_debye[1] << " Debye\n";
    std::cout << "  mu_z = " << mu_debye[2] << " Debye\n";
    std::cout << "  |mu| = " << magnitude_debye << " Debye\n";
    
    std::cout << "------------------------------------------\n";
    std::cout << "Reference values (experimental):\n";
    std::cout << "  HCl:  1.08 Debye\n";
    std::cout << "  H2O:  1.85 Debye\n";
    std::cout << "  H2:   0.00 Debye\n";
}

void DipoleMoment::write(const std::string& filename) const {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("DipoleMoment: could not write to: " + filename);
    }

    file << "# Dipole moment in Debye: mu_x mu_y mu_z magnitude\n";
    file << mu_debye[0] << " "
         << mu_debye[1] << " "
         << mu_debye[2] << " "
         << magnitude_debye << "\n";
}

std::array<double, 3> DipoleMoment::get_components_debye() const {
    return mu_debye;
}

double DipoleMoment::get_magnitude_debye() const {
    return magnitude_debye;
}
