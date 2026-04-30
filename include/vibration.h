#ifndef VIBRATION_H
#define VIBRATION_H

#include "hessian.h"
#include "molecule.h"

#include <string>
#include <vector>

class Vibrations {
private:
    std::vector<double> frequencies;

    bool is_linear_molecule(const Molecule& molecule) const;

public:
    void compute(const Molecule& molecule, const Hessian& hessian);

    void print_frequencies() const;
    void write_frequencies(const std::string& filename) const;

    std::vector<double> get_frequencies() const;
};

#endif
