//
#include "molecule.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

/*
    Simple CNDO-energy executable placeholder.

    Purpose:
    - takes an xyz file
    - computes a simple geometry-dependent energy
    - prints one number

    This lets the finite-difference Hessian pipeline run.

    Later, this file can be replaced with the real CNDO/2 SCF energy code.
*/

double get_atomic_number(const std::string& symbol) {
    if (symbol == "H") return 1.0;
    if (symbol == "C") return 6.0;
    if (symbol == "N") return 7.0;
    if (symbol == "O") return 8.0;
    if (symbol == "F") return 9.0;
    if (symbol == "Cl") return 17.0;

    return 1.0;
}

double compute_energy(const Molecule& molecule) {
    std::vector<std::string> symbols = molecule.get_symbols();
    std::vector<double> coords = molecule.get_coordinates();

    int n = molecule.get_num_atoms();

    double energy = 0.0;

    for (int i = 0; i < n; i++) {
        double xi = coords[3 * i];
        double yi = coords[3 * i + 1];
        double zi = coords[3 * i + 2];

        double zi_charge = get_atomic_number(symbols[i]);

        for (int j = i + 1; j < n; j++) {
            double xj = coords[3 * j];
            double yj = coords[3 * j + 1];
            double zj = coords[3 * j + 2];

            double zj_charge = get_atomic_number(symbols[j]);

            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;

            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r > 1.0e-8) {
                // Simple nuclear-repulsion-like term.
                energy += (zi_charge * zj_charge) / r;
            }
        }
    }

    return energy;
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cout << "Usage: ./cndo_energy molecule.xyz\n";
            return 1;
        }

        Molecule molecule;
        molecule.read_xyz(argv[1]);

        double energy = compute_energy(molecule);

        // CNDOEngine reads the last number printed, so keep this simple.
        std::cout << energy << std::endl;
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
