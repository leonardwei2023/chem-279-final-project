#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>

class Molecule {
private:
    int num_atoms;
    std::vector<std::string> symbols;
    // x1 y1 z1 x2 y2 z2 ...
    std::vector<double> coordinates; 
    std::vector<double> masses;

    double get_mass_from_symbol(const std::string& symbol) const;

public:
    Molecule();

    void read_xyz(const std::string& filename);

    int get_num_atoms() const;
    std::vector<std::string> get_symbols() const;
    std::vector<double> get_coordinates() const;
    std::vector<double> get_masses() const;
};

#endif
