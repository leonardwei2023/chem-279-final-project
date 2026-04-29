#ifndef VIBRATIONAL_FREQUENCY_H
#define VIBRATIONAL_FREQUENCY_H

#include <armadillo>
#include "molecule.h"

class VibrationalFrequencyAnalyzer {
public:
    VibrationalFrequencyAnalyzer(double step_size);

    arma::mat compute_hessian(Molecule mol);
    arma::mat mass_weight_hessian(const arma::mat& hessian, const Molecule& mol);
    arma::vec compute_frequencies(const arma::mat& mass_weighted_hessian);
    void print_frequencies(const arma::vec& frequencies);

private:
    double h;
};

#endif
