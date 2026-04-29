#include <armadillo>
#include <iostream>

#include "molecule.h"
#include "vibrational_frequency.h"

int main()
{
    std::cout << "Starting vibrational frequency test...\n";

    Molecule mol;
    mol.set_num_atoms(2);

    // Simple H2 test molecule.
    // Coordinates are in Bohr.
    
    mol.set_symbol(0, "H");
    mol.set_symbol(1, "H");

    mol.set_atomic_mass(0, 1.00784);
    mol.set_atomic_mass(1, 1.00784);

    mol.coordinates(0, 0) = 0.0;
    mol.coordinates(0, 1) = 0.0;
    mol.coordinates(0, 2) = 0.0;

    mol.coordinates(1, 0) = 0.0;
    mol.coordinates(1, 1) = 0.0;
    mol.coordinates(1, 2) = 1.4;

    double step_size = 0.005;

    VibrationalFrequencyAnalyzer vib(step_size);

    arma::mat hessian = vib.compute_hessian(mol);
    arma::mat mass_weighted_hessian = vib.mass_weight_hessian(hessian, mol);
    arma::vec frequencies = vib.compute_frequencies(mass_weighted_hessian);

    vib.print_frequencies(frequencies);

    std::cout << "\nDone.\n";

    return 0;
}
