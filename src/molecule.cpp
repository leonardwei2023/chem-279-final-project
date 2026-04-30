#include "molecule.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

Molecule::Molecule() {
    num_atoms = 0;
}

double Molecule::get_mass_from_symbol(const std::string& symbol) const {
    if (symbol == "H") return 1.0079;
    if (symbol == "C") return 12.011;
    if (symbol == "N") return 14.007;
    if (symbol == "O") return 15.999;
    if (symbol == "F") return 18.998;
    if (symbol == "Cl") return 35.45;

    throw std::runtime_error("Unknown atom symbol: " + symbol);
}

void Molecule::read_xyz(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open xyz file: " + filename);
    }

    file >> num_atoms;

    std::string line;
    std::getline(file, line); 
    std::getline(file, line); 

    symbols.clear();
    coordinates.clear();
    masses.clear();

    for (int i = 0; i < num_atoms; i++) {
        std::string symbol;
        double x, y, z;

        file >> symbol >> x >> y >> z;

        symbols.push_back(symbol);
        coordinates.push_back(x);
        coordinates.push_back(y);
        coordinates.push_back(z);
        masses.push_back(get_mass_from_symbol(symbol));
    }

    file.close();
}

int Molecule::get_num_atoms() const {
    return num_atoms;
}

std::vector<std::string> Molecule::get_symbols() const {
    return symbols;
}

std::vector<double> Molecule::get_coordinates() const {
    return coordinates;
}

std::vector<double> Molecule::get_masses() const {
    return masses;
}
