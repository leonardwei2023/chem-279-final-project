#include "molecule.h"
#include "hessian.h"
#include "vibration.h"

// later: #include "dipole.h"

#include <iostream>
#include <string>
#include <exception>


int main(int argc, char* argv[]) {
    try {
        if (argc < 3) {
            std::cout << "Usage:\n";
            std::cout << "./program vibration xyz hessian\n";
            std::cout << "./program dipole xyz\n";
            return 1;
        }

        std::string mode = argv[1];

        if (mode == "vibration") {
            if (argc != 4) {
                std::cout << "Usage: ./program vibration xyz hessian\n";
                return 1;
            }

            std::string xyz_file = argv[2];
            std::string hessian_file = argv[3];

            Molecule molecule;
            molecule.read_xyz(xyz_file);

            int size = 3 * molecule.get_num_atoms();

            Hessian hessian;
            hessian.read_from_file(hessian_file, size);

            Vibrations vibrations;
            vibrations.compute(molecule, hessian);
            vibrations.print_frequencies();
        }
        else if (mode == "dipole") {
            std::string xyz_file = argv[2];

            Molecule molecule;
            molecule.read_xyz(xyz_file);

            // placeholder for now

            
            std::cout << "Dipole calculation not implemented yet.\n";

            // later: Add these for dipole moment:
            // Dipole dipole;
            // dipole.compute(molecule);
            // dipole.print();
        }
        else {
            std::cout << "Unknown mode: " << mode << std::endl;
        }
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
