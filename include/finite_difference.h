#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include <armadillo>
#include "molecule.h"

// Builds the Cartesian Hessian using finite differences of the CNDO/2 energy:

arma::mat compute_hessian_finite_difference(Molecule mol, double step_size);

// Converts the normal Hessian into a mass-weighted Hessian:

arma::mat mass_weight_hessian(const arma::mat& hessian, const Molecule& mol);

// Diagonalizes the mass-weighted Hessian and converts eigenvalues to cm^-1:

arma::vec compute_frequencies_cm(const arma::mat& mass_weighted_hessian);

// Prints all modes in a clean table:

void print_frequencies(const arma::vec& frequencies);

#endif
