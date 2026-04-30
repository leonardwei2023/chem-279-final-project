
#ifndef VIBRATIONS_H
#define VIBRATIONS_H

#include "molecule.h"
#include "hessian.h"

#include <vector>

class Vibrations {
private:
    std::vector<double> frequencies;

public:
    void compute(const Molecule& molecule, const Hessian& hessian);
    void print_frequencies() const;
};

#endif
