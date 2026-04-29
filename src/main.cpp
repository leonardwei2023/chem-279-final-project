
#include "vibrational_frequency.h"



double step_size = 0.005;

VibrationalFrequencyAnalyzer vib(step_size);

arma::mat hessian = vib.compute_hessian(mol);
arma::mat mass_weighted_hessian = vib.mass_weight_hessian(hessian, mol);
arma::vec frequencies = vib.compute_frequencies(mass_weighted_hessian);

vib.print_frequencies(frequencies);
