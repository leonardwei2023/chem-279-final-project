#ifndef VIBRATIONAL_FREQUENCY_H
#define VIBRATIONAL_FREQUENCY_H

#include <armadillo>
#include "molecule.h"

// Main class for finite-difference vibrational frequency analysis:

class VibrationalFrequencyAnalyzer {
public:

    // Constructor:

    VibrationalFrequencyAnalyzer(double step_size);

    // Builds the Cartesian Hessian using finite differences:

    arma::mat compute_hessian(Molecule mol);

    // Converts Hessian to mass-weighted Hessian:

    arma::mat mass_weight_hessian(const arma::mat& hessian, const Molecule& mol);

    // Diagonalizes mass-weighted Hessian and returns frequencies in cm^-1:

    arma::vec compute_frequencies(const arma::mat& mass_weighted_hessian);

    // Prints frequencies nicely:

    void print_frequencies(const arma::vec& frequencies);

private:
    double h;
};

#endif
