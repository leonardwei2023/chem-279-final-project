#include "molecule.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

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
        double Zi = get_atomic_number(symbols[i]);

        for (int j = i + 1; j < n; j++) {
            double xj = coords[3 * j];
            double yj = coords[3 * j + 1];
            double zj = coords[3 * j + 2];
            double Zj = get_atomic_number(symbols[j]);

            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;

            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r > 1.0e-8) {
                energy += (Zi * Zj) / r;
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

        std::cout << energy << std::endl;
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
