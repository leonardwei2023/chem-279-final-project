#include "cndo_engine.h"
#include "dipole.h"
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

static std::string get_molecule_name(const std::string& filepath) {
    size_t slash = filepath.find_last_of("/\\");
    std::string filename = (slash == std::string::npos)
                           ? filepath
                           : filepath.substr(slash + 1);
    size_t dot = filename.find_last_of(".");
    if (dot != std::string::npos) filename = filename.substr(0, dot);
    return filename;
}

static void print_molecule_info(const std::string& xyz_file,
                                const Molecule& molecule) {
    std::cout << "\n============================\n";
    std::cout << "Molecule : " << get_molecule_name(xyz_file) << "\n";
    std::cout << "File     : " << xyz_file << "\n";
    std::cout << "Atoms    : " << molecule.get_num_atoms() << "\n";
    std::cout << "============================\n";
}

static int get_expected_modes(const Molecule& molecule) {
    int n = molecule.get_num_atoms();
    return (n == 2) ? 3 * n - 5 : 3 * n - 6;
}

static void print_usage() {
    std::cout << "\nUsage:\n"
              << "  ./vibrational_frequency vibration        <xyz> <hessian_dat> [--animate]\n"
              << "  ./vibrational_frequency finite-diff      <xyz> <out_hessian> <step>\n"
              << "  ./vibrational_frequency finite-diff-vib  <xyz> <out_hessian> <step> [--animate]\n"
              << "  ./vibrational_frequency dipole           <xyz>\n"
              << "  ./vibrational_frequency dipole-scf       <xyz> <p_diagonal_file>\n"
              << "  ./vibrational_frequency validate         <computed_freq> <reference_freq>\n"
              << "\n"
              << "Environment variable for finite-diff modes:\n"
              << "  export CNDO_ENERGY_CMD=\"./cndo_energy\"\n"
              << "\n"
              << "Options:\n"
              << "  --animate   also write normal_modes.xyz for Avogadro/VMD visualization\n";
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        std::string mode = argv[1];

        if (mode == "vibration") {
            if (argc < 4) { print_usage(); return 1; }

            std::string xyz_file = argv[2];
            std::string hessian_file = argv[3];
            bool animate = (argc >= 5 && std::string(argv[4]) == "--animate");

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

            if (animate) {
                vibrations.write_normal_mode_xyz(molecule, "normal_modes.xyz");
            }
        }

        else if (mode == "finite-diff") {
            if (argc != 5) { print_usage(); return 1; }

            std::string xyz_file = argv[2];
            std::string output_hessian = argv[3];
            double step = std::stod(argv[4]);

            const char* cmd = std::getenv("CNDO_ENERGY_CMD");
            if (!cmd) {
                std::cout << "Error: CNDO_ENERGY_CMD not set.\n"
                          << "Example: export CNDO_ENERGY_CMD=\"./cndo_energy\"\n";
                return 1;
            }

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            CNDOEngine engine(cmd);
            FiniteDifference fd(step);

            Hessian hessian = fd.compute_hessian(molecule, engine);
            hessian.write_to_file(output_hessian);

            std::cout << "Hessian written to: " << output_hessian << "\n";
        }

        else if (mode == "finite-diff-vib") {
            if (argc < 5) { print_usage(); return 1; }

            std::string xyz_file = argv[2];
            std::string output_hessian = argv[3];
            double step = std::stod(argv[4]);
            bool animate = (argc >= 6 && std::string(argv[5]) == "--animate");

            const char* cmd = std::getenv("CNDO_ENERGY_CMD");
            if (!cmd) {
                std::cout << "Error: CNDO_ENERGY_CMD not set.\n"
                          << "Example: export CNDO_ENERGY_CMD=\"./cndo_energy\"\n";
                return 1;
            }

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            CNDOEngine engine(cmd);
            FiniteDifference fd(step);

            Hessian hessian = fd.compute_hessian(molecule, engine);
            hessian.write_to_file(output_hessian);

            std::cout << "Hessian written to: " << output_hessian << "\n";
            std::cout << "Expected vibrational modes: "
                      << get_expected_modes(molecule) << "\n";

            Vibrations vibrations;
            vibrations.compute(molecule, hessian);
            vibrations.print_frequencies();
            vibrations.write_frequencies("computed_frequencies.dat");

            if (animate) {
                vibrations.write_normal_mode_xyz(molecule, "normal_modes.xyz");
            }
        }

        else if (mode == "dipole") {
            if (argc != 3) { print_usage(); return 1; }

            std::string xyz_file = argv[2];

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            std::cout << "Running CNDO/2 SCF to generate p_diagonal.dat...\n";

            CNDOEngine engine("./cndo_energy");
            double energy = engine.compute_energy(molecule);

            std::cout << "SCF energy = " << energy << " Hartree\n";

            std::vector<double> p_diagonal =
                Validation::read_frequencies("p_diagonal.dat");

            std::cout << "Loaded " << p_diagonal.size()
                      << " SCF atomic populations from p_diagonal.dat\n";

            DipoleMoment dipole;
            dipole.compute(molecule, p_diagonal);
            dipole.print();
            dipole.write("dipole_moment.dat");
        }

        else if (mode == "dipole-scf") {
            if (argc != 4) { print_usage(); return 1; }

            std::string xyz_file = argv[2];
            std::string pdiag_file = argv[3];

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            std::vector<double> p_diagonal =
                Validation::read_frequencies(pdiag_file);

            std::cout << "Loaded " << p_diagonal.size()
                      << " diagonal density values from: " << pdiag_file << "\n";

            DipoleMoment dipole;
            dipole.compute(molecule, p_diagonal);
            dipole.print();
            dipole.write("dipole_moment.dat");
        }

        else if (mode == "validate") {
            if (argc != 4) { print_usage(); return 1; }

            std::vector<double> computed =
                Validation::read_frequencies(argv[2]);
            std::vector<double> reference =
                Validation::read_frequencies(argv[3]);

            Validation::compare(computed, reference);
        }

        else {
            std::cout << "Unknown mode: " << mode << "\n";
            print_usage();
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
