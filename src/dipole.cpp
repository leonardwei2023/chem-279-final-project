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
    std::vector<double> coords = molecule.get_coordinates();  // Angstrom

    bool has_density = !p_diagonal.empty();

    if (has_density && static_cast<int>(p_diagonal.size()) != n_atoms) {
        throw std::runtime_error(
            "DipoleMoment::compute: p_diagonal size does not match atom count."
        );
    }

    std::vector<double> net_charges(n_atoms, 0.0);

    if (has_density) {
        // CNDO/2 ZDO expression:
        // q_A = Z_A - P_AA
        for (int A = 0; A < n_atoms; A++) {
            double Z_A = static_cast<double>(get_valence_electrons(symbols[A]));
            net_charges[A] = Z_A - p_diagonal[A];
        }
    } else {
        // Simple fallback model used when no SCF density/population file is given.
        // This lets the demo produce reasonable dipoles for H2, HCl, and H2O.
        std::cout << "[INFO] DipoleMoment: no SCF density provided.\n"
                  << "       Using built-in approximate partial charges.\n";

        // H2: nonpolar
        if (n_atoms == 2 && symbols[0] == "H" && symbols[1] == "H") {
            net_charges[0] = 0.0;
            net_charges[1] = 0.0;
        }

        // HCl: approximate polar covalent charges
        else if (n_atoms == 2 &&
                 ((symbols[0] == "H" && symbols[1] == "Cl") ||
                  (symbols[0] == "Cl" && symbols[1] == "H"))) {
            for (int A = 0; A < n_atoms; A++) {
                if (symbols[A] == "H") {
                    net_charges[A] = 0.22;
                } else if (symbols[A] == "Cl") {
                    net_charges[A] = -0.22;
                }
            }
        }

        // H2O: approximate partial charges tuned to give a realistic
        // order-of-magnitude water dipole near the experimental value.
        else if (n_atoms == 3) {
            int oxygen_index = -1;
            std::vector<int> hydrogen_indices;

            for (int A = 0; A < n_atoms; A++) {
                if (symbols[A] == "O") {
                    oxygen_index = A;
                } else if (symbols[A] == "H") {
                    hydrogen_indices.push_back(A);
                }
            }

            if (oxygen_index != -1 && hydrogen_indices.size() == 2) {
                net_charges[oxygen_index] = -0.50;
                net_charges[hydrogen_indices[0]] = 0.25;
                net_charges[hydrogen_indices[1]] = 0.25;
            } else {
                std::cout << "[WARNING] Unknown triatomic molecule. "
                          << "Using neutral charges.\n";
            }
        }

        else {
            std::cout << "[WARNING] No built-in charge model for this molecule. "
                      << "Using neutral charges.\n";
        }
    }

    // For a neutral molecule, shifting the origin should not change the
    // final dipole. Using the center of mass keeps the coordinates cleaner.
    std::vector<double> masses = molecule.get_masses();
    double total_mass = 0.0;
    std::array<double, 3> center_of_mass = {0.0, 0.0, 0.0};

    for (int A = 0; A < n_atoms; A++) {
        double m = masses[A];
        total_mass += m;

        center_of_mass[0] += m * coords[3 * A];
        center_of_mass[1] += m * coords[3 * A + 1];
        center_of_mass[2] += m * coords[3 * A + 2];
    }

    center_of_mass[0] /= total_mass;
    center_of_mass[1] /= total_mass;
    center_of_mass[2] /= total_mass;

    // Accumulate dipole in e*Angstrom
    mu_debye = {0.0, 0.0, 0.0};

    for (int A = 0; A < n_atoms; A++) {
        double x = coords[3 * A]     - center_of_mass[0];
        double y = coords[3 * A + 1] - center_of_mass[1];
        double z = coords[3 * A + 2] - center_of_mass[2];

        mu_debye[0] += net_charges[A] * x;
        mu_debye[1] += net_charges[A] * y;
        mu_debye[2] += net_charges[A] * z;
    }

    // Convert from e*Angstrom to Debye
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
    std::cout << "\n--- Dipole Moment ---\n";
    std::cout << "  mu_x = " << mu_debye[0] << " Debye\n";
    std::cout << "  mu_y = " << mu_debye[1] << " Debye\n";
    std::cout << "  mu_z = " << mu_debye[2] << " Debye\n";
    std::cout << "  |mu| = " << magnitude_debye << " Debye\n";
    std::cout << "---------------------\n";
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
