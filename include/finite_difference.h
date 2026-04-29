
#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include <armadillo>
#include "molecule.h"

arma::mat compute_hessian_finite_difference(Molecule mol, double h);
arma::mat mass_weight_hessian(const arma::mat& H, const Molecule& mol);
arma::vec compute_frequencies_cm(const arma::mat& Hmw);

void print_frequencies(const arma::vec& freqs);

#endif
