
#include <armadillo>
#include <iostream>

#include "molecule.h"
#include "vibrational_frequency.h"

int main()
{
    std::cout << "Starting Vibrational Frequency Analysis...\n";

    // Create a simple test molecule (H2)
    // This lets the code run even without file input.

  
    Molecule mol;

    mol.set_num_atoms(2);

    // Coordinates in Bohr 
  
    mol.coordinates = arma::mat(2, 3, arma::fill::zeros);

    // Place two hydrogen atoms along z-axis
  
    mol.coordinates(0, 2) = 0.0;
    mol.coordinates(1, 2) = 1.4;  // bond length (approx)

    // Atomic masses (amu)
  
    mol.set_atomic_mass(0, 1.0);  // H
    mol.set_atomic_mass(1, 1.0);  // H


      // Vibrational frequency analysis

      double step_size = 0.005;

    VibrationalFrequencyAnalyzer vib(step_size);

    arma::mat hessian = vib.compute_hessian(mol);
    arma::mat mass_weighted_hessian = vib.mass_weight_hessian(hessian, mol);
    arma::vec frequencies = vib.compute_frequencies(mass_weighted_hessian);

    vib.print_frequencies(frequencies);

    std::cout << "\nDone.\n";

    return 0;
}
