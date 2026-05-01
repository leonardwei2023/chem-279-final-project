std::cout << "NEW ENERGY MODEL RUNNING\n";
#include "cndo_engine.h"
#include "finite_difference.h"
#include "hessian.h"
#include "molecule.h"
#include "validation.h"
#include "vibration.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

std::string get_molecule_name(const std::string& filepath) {
    size_t slash_pos = filepath.find_last_of("/\\");
    std::string filename;

    if (slash_pos == std::string::npos) {
        filename = filepath;
    }
    else {
        filename = filepath.substr(slash_pos + 1);
    }

    size_t dot_pos = filename.find_last_of(".");
    if (dot_pos != std::string::npos) {
        filename = filename.substr(0, dot_pos);
    }

    return filename;
}

void print_molecule_info(const std::string& xyz_file, const Molecule& molecule) {
    std::cout << "\n============================\n";
    std::cout << "Molecule: " << get_molecule_name(xyz_file) << "\n";
    std::cout << "File: " << xyz_file << "\n";
    std::cout << "Atoms: " << molecule.get_num_atoms() << "\n";
    std::cout << "============================\n";
}

int get_expected_modes(const Molecule& molecule) {
    int n = molecule.get_num_atoms();

    if (n == 2) {
        return 3 * n - 5;
    }

    return 3 * n - 6;
}

void print_usage() {
    std::cout << "Usage:\n";
    std::cout << "./vibrational_frequency vibration xyz_file hessian_file\n";
    std::cout << "./vibrational_frequency finite-diff xyz_file output_hessian step\n";
    std::cout << "./vibrational_frequency finite-diff-vib xyz_file output_hessian step\n";
    std::cout << "./vibrational_frequency validate computed_freq reference_freq\n";
    std::cout << "./vibrational_frequency dipole xyz_file\n";
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        std::string mode = argv[1];

        if (mode == "vibration") {
            if (argc != 4) {
                print_usage();
                return 1;
            }

            std::string xyz_file = argv[2];
            std::string hessian_file = argv[3];

            Molecule molecule;
            molecule.read_xyz(xyz_file);

            print_molecule_info(xyz_file, molecule);

            Hessian hessian;
            hessian.read_from_file(hessian_file, molecule.get_num_coordinates());

            std::cout << "Expected vibrational modes: "
                      << get_expected_modes(molecule) << "\n";

            Vibrations vibrations;
            vibrations.compute(molecule, hessian);
            vibrations.print_frequencies();
            vibrations.write_frequencies("computed_frequencies.dat");
        }

        else if (mode == "finite-diff") {
            if (argc != 5) {
                print_usage();
                return 1;
            }

            std::string xyz_file = argv[2];
            std::string output_hessian = argv[3];
            double step = std::stod(argv[4]);

            const char* command_from_env = std::getenv("CNDO_ENERGY_CMD");

            if (command_from_env == nullptr) {
                std::cout << "Error: CNDO_ENERGY_CMD is not set.\n";
                std::cout << "Example:\n";
                std::cout << "export CNDO_ENERGY_CMD=\"./build/cndo_energy\"\n";
                return 1;
            }

            Molecule molecule;
            molecule.read_xyz(xyz_file);

            print_molecule_info(xyz_file, molecule);

            CNDOEngine engine(command_from_env);
            FiniteDifference finite_difference(step);

            Hessian hessian = finite_difference.compute_hessian(molecule, engine);
            hessian.write_to_file(output_hessian);

            std::cout << "Hessian written to: " << output_hessian << "\n";
        }

        else if (mode == "finite-diff-vib") {
            if (argc != 5) {
                print_usage();
                return 1;
            }

            std::string xyz_file = argv[2];
            std::string output_hessian = argv[3];
            double step = std::stod(argv[4]);

            const char* command_from_env = std::getenv("CNDO_ENERGY_CMD");

            if (command_from_env == nullptr) {
                std::cout << "Error: CNDO_ENERGY_CMD is not set.\n";
                std::cout << "Example:\n";
                std::cout << "export CNDO_ENERGY_CMD=\"./build/cndo_energy\"\n";
                return 1;
            }

            Molecule molecule;
            molecule.read_xyz(xyz_file);

            print_molecule_info(xyz_file, molecule);

            CNDOEngine engine(command_from_env);
            FiniteDifference finite_difference(step);

            Hessian hessian = finite_difference.compute_hessian(molecule, engine);
            hessian.write_to_file(output_hessian);

            std::cout << "Hessian written to: " << output_hessian << "\n";
            std::cout << "Expected vibrational modes: "
                      << get_expected_modes(molecule) << "\n";

            Vibrations vibrations;
            vibrations.compute(molecule, hessian);
            vibrations.print_frequencies();
            vibrations.write_frequencies("computed_frequencies.dat");
        }

        else if (mode == "validate") {
            if (argc != 4) {
                print_usage();
                return 1;
            }

            std::vector<double> computed = Validation::read_frequencies(argv[2]);
            std::vector<double> reference = Validation::read_frequencies(argv[3]);

            Validation::compare(computed, reference);
        }

        else if (mode == "dipole") {
            if (argc != 3) {
                print_usage();
                return 1;
            }

            Molecule molecule;
            molecule.read_xyz(argv[2]);

            print_molecule_info(argv[2], molecule);

            std::cout << "Dipole moment calculation placeholder.\n";
            std::cout << "Dipole moment is not implemented yet.\n";
        }

        else {
            std::cout << "Unknown mode: " << mode << "\n";
            print_usage();
            return 1;
        }
    }
    catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
