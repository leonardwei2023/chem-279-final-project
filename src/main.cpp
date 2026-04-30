#include "molecule.h"
#include "hessian.h"
#include "vibration.h"

#include <iostream>

int main() {
    try {
        Molecule molecule;
        molecule.read_xyz("input/h2o.xyz");

        int hessian_size = 3 * molecule.get_num_atoms();

        Hessian hessian;
        hessian.read_from_file("input/hessian.dat", hessian_size);

        Vibrations vibrations;
        vibrations.compute(molecule, hessian);
        vibrations.print_frequencies();
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
