#include "dipole.h"
#include "molecule.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

// -----------------------------------------------------------------------
// CNDO/2 valence electron counts (only valence electrons are treated)
// These match the CNDO/2 parameterization from Lecture 13.
//   H:  1 valence electron  (1s)
//   C:  4 valence electrons (2s, 2px, 2py, 2pz)
//   N:  5 valence electrons
//   O:  6 valence electrons
//   F:  7 valence electrons
//   Cl: 7 valence electrons (3s, 3p)
// -----------------------------------------------------------------------
int DipoleMoment::get_valence_electrons(const std::string& symbol) {
    if (symbol == "H")  return 1;
    if (symbol == "C")  return 4;
    if (symbol == "N")  return 5;
    if (symbol == "O")  return 6;
    if (symbol == "F")  return 7;
    if (symbol == "Cl") return 7;
    throw std::runtime_error("DipoleMoment: unknown element symbol: " + symbol);
}

// -----------------------------------------------------------------------
// Compute dipole moment
//
// Formula (Lecture 16, slide 10):
//   d_SCF = sum_{mu,nu} P_{mu,nu} * d_{mu,nu}
//
// Under CNDO/2 ZDO approximation, the only surviving dipole integrals
// are on-site:  d_{mu,mu} = R_A  (position of atom A owning orbital mu)
//
// So the full expression becomes:
//   mu_j = sum_A [ Z_A - (sum_{mu on A} P_{mu,mu}) ] * R_A^j
//
// where j in {x, y, z}, Z_A is the CNDO/2 valence charge, and
// sum_{mu on A} P_{mu,mu} is the electron population on atom A.
//
// If p_diagonal is empty (no SCF available yet), we fall back to
// a pure nuclear point-charge model as a placeholder, printing a warning.
// -----------------------------------------------------------------------
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
            "DipoleMoment::compute: p_diagonal size (" +
            std::to_string(p_diagonal.size()) +
            ") does not match number of atoms (" +
            std::to_string(n_atoms) + ")."
        );
    }

    if (!has_density) {
        std::cout << "[WARNING] DipoleMoment: no density matrix provided.\n"
                  << "          Using neutral-atom approximation (all net charges = 0).\n"
                  << "          Connect CNDO/2 SCF output for proper results.\n";
    }

    // Accumulate dipole vector in e*Angstrom
    mu_debye = {0.0, 0.0, 0.0};

    for (int A = 0; A < n_atoms; A++) {
        double Z_A = static_cast<double>(get_valence_electrons(symbols[A]));

        // Electron population on atom A from diagonal density
        // In neutral-atom fallback, assume electrons exactly cancel nuclear charge
        double elec_pop_A = has_density ? p_diagonal[A] : Z_A;

        // Net charge on atom A (in units of e)
        double net_charge_A = Z_A - elec_pop_A;

        double x_A = coords[3 * A];
        double y_A = coords[3 * A + 1];
        double z_A = coords[3 * A + 2];

        mu_debye[0] += net_charge_A * x_A;
        mu_debye[1] += net_charge_A * y_A;
        mu_debye[2] += net_charge_A * z_A;
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
    std::cout << "\n--- Dipole Moment (CNDO/2 ZDO approximation) ---\n";
    std::cout << "  mu_x = " << mu_debye[0] << " Debye\n";
    std::cout << "  mu_y = " << mu_debye[1] << " Debye\n";
    std::cout << "  mu_z = " << mu_debye[2] << " Debye\n";
    std::cout << "  |mu| = " << magnitude_debye << " Debye\n";
    std::cout << "------------------------------------------------\n";
    std::cout << "Reference values (experimental):\n";
    std::cout << "  HCl:  1.08 Debye\n";
    std::cout << "  H2O:  1.85 Debye\n";
    std::cout << "  H2:   0.00 Debye (nonpolar, expected)\n";
}

void DipoleMoment::write(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("DipoleMoment: could not write to: " + filename);
    }

    file << "# Dipole moment (Debye): mu_x  mu_y  mu_z  |mu|\n";
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
