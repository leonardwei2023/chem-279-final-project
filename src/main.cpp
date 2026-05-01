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

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

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
              << "  export CNDO_ENERGY_CMD=\"./build/cndo_energy\"\n"
              << "\n"
              << "Options:\n"
              << "  --animate   also write normal_modes.xyz for Avogadro/VMD visualization\n";
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------
int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        std::string mode = argv[1];

        // ----------------------------------------------------------------
        // MODE: vibration
        //   Load a pre-computed Hessian and compute frequencies.
        //   ./vibrational_frequency vibration h2.xyz h2_hessian.dat [--animate]
        // ----------------------------------------------------------------
        if (mode == "vibration") {
            if (argc < 4) { print_usage(); return 1; }

            std::string xyz_file     = argv[2];
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
                vibrations.write_normal_mode_xyz(
                    molecule, "normal_modes.xyz"
                );
            }
        }

        // ----------------------------------------------------------------
        // MODE: finite-diff
        //   Compute Hessian via finite differences using CNDO_ENERGY_CMD.
        // ----------------------------------------------------------------
        else if (mode == "finite-diff") {
            if (argc != 5) { print_usage(); return 1; }

            std::string xyz_file      = argv[2];
            std::string output_hessian = argv[3];
            double step               = std::stod(argv[4]);

            const char* cmd = std::getenv("CNDO_ENERGY_CMD");
            if (!cmd) {
                std::cout << "Error: CNDO_ENERGY_CMD not set.\n"
                          << "Example: export CNDO_ENERGY_CMD=\"./build/cndo_energy\"\n";
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

        // ----------------------------------------------------------------
        // MODE: finite-diff-vib
        //   Compute Hessian AND then compute/print frequencies.
        // ----------------------------------------------------------------
        else if (mode == "finite-diff-vib") {
            if (argc < 5) { print_usage(); return 1; }

            std::string xyz_file      = argv[2];
            std::string output_hessian = argv[3];
            double step               = std::stod(argv[4]);
            bool animate = (argc >= 6 && std::string(argv[5]) == "--animate");

            const char* cmd = std::getenv("CNDO_ENERGY_CMD");
            if (!cmd) {
                std::cout << "Error: CNDO_ENERGY_CMD not set.\n"
                          << "Example: export CNDO_ENERGY_CMD=\"./build/cndo_energy\"\n";
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
                vibrations.write_normal_mode_xyz(
                    molecule, "normal_modes.xyz"
                );
            }
        }

        // ----------------------------------------------------------------
        // MODE: dipole
        //   Compute dipole moment using CNDO/2 ZDO formula.
        //   Without SCF density, uses neutral-atom fallback (mu ~ 0).
        //   ./vibrational_frequency dipole h2.xyz
        // ----------------------------------------------------------------
        else if (mode == "dipole") {
            if (argc != 3) { print_usage(); return 1; }

            std::string xyz_file = argv[2];

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            // No density matrix available yet — use neutral atom fallback.
            // This will give mu = 0 for all molecules since net charges = 0.
            // Once your CNDO/2 SCF is connected, use 'dipole-scf' mode below.
            std::vector<double> p_diagonal;  // empty = neutral atom fallback

            DipoleMoment dipole;
            dipole.compute(molecule, p_diagonal);
            dipole.print();
            dipole.write("dipole_moment.dat");
        }

        // ----------------------------------------------------------------
        // MODE: dipole-scf
        //   Compute dipole moment using diagonal of CNDO/2 density matrix.
        //
        //   Formula (Lecture 16, slide 10):
        //     mu_j = sum_A (Z_A - P_AA) * R_A^j   [e*Ang] * EA_TO_DEBYE
        //
        //   The p_diagonal file has one value per atom:
        //     p_diagonal[A] = sum_{mu on atom A} P_{mu,mu}
        //   i.e. the electron population on each atom from your SCF code.
        //
        //   ./vibrational_frequency dipole-scf h2.xyz p_diagonal.dat
        // ----------------------------------------------------------------
        else if (mode == "dipole-scf") {
            if (argc != 4) { print_usage(); return 1; }

            std::string xyz_file    = argv[2];
            std::string pdiag_file  = argv[3];

            Molecule molecule;
            molecule.read_xyz(xyz_file);
            print_molecule_info(xyz_file, molecule);

            // Read per-atom electron populations from file
            std::vector<double> p_diagonal =
                Validation::read_frequencies(pdiag_file);  // reuses file reader

            std::cout << "Loaded " << p_diagonal.size()
                      << " diagonal density values from: " << pdiag_file << "\n";

            DipoleMoment dipole;
            dipole.compute(molecule, p_diagonal);
            dipole.print();
            dipole.write("dipole_moment.dat");
        }

        // ----------------------------------------------------------------
        // MODE: validate
        //   Compare computed vs reference frequencies.
        // ----------------------------------------------------------------
        else if (mode == "validate") {
            if (argc != 4) { print_usage(); return 1; }

            std::vector<double> computed  = Validation::read_frequencies(argv[2]);
            std::vector<double> reference = Validation::read_frequencies(argv[3]);

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
