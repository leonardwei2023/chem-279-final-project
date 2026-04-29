#include "molecule.h"

Molecule::Molecule()
{
}

void Molecule::set_num_atoms(int n)
{
    coordinates = arma::mat(n, 3, arma::fill::zeros);
    symbols.resize(n);
    masses.resize(n, 1.0);
}

int Molecule::num_atoms() const
{
    return coordinates.n_rows;
}

void Molecule::set_symbol(int index, const std::string& symbol)
{
    symbols[index] = symbol;
}

std::string Molecule::atomic_symbol(int index) const
{
    return symbols[index];
}

void Molecule::set_atomic_mass(int index, double mass)
{
    masses[index] = mass;
}

double Molecule::atomic_mass(int index) const
{
    return masses[index];
}
