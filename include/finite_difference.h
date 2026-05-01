
#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include "cndo_engine.h"
#include "hessian.h"
#include "molecule.h"

class FiniteDifference {
private:
    double step_size;

public:
    FiniteDifference(double step);

    Hessian compute_hessian(const Molecule& molecule, CNDOEngine& engine);
};

#endif
