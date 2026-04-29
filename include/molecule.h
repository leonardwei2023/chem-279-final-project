
#ifndef MOLECULE_H
#define MOLECULE_H

#include <armadillo>
#include <vector>
#include <string>

class Molecule {
public:
    arma::mat coordinates;
    std::vector<std::string> symbols;
    std::vector<double> masses;

    Molecule();

    void set_num_atoms(int n);
    int num_atoms() const;

    void set_symbol(int index, const std::string& symbol);
    std::string atomic_symbol(int index) const;

    void set_atomic_mass(int index, double mass);
    double atomic_mass(int index) const;
};

#endif
