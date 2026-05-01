#include "vibration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

// Hessian from finite_difference.cpp should be in Hartree / Bohr^2.
// For mass-weighting in atomic units, convert masses:
// 1 amu = 1822.888486 electron masses.


static constexpr double AMU_TO_AU = 1822.888486;

// Converts sqrt(eigenvalue in atomic units) to cm^-1:

static constexpr double AU_TO_CM = 5140.48;

// Small frequencies are treated as translation/rotation modes:

static constexpr double FREQ_ZERO_THRESHOLD = 10.0;

bool Vibrations::is_linear_molecule(const Molecule& molecule) const {
    return molecule.get_num_atoms() == 2;
}

void Vibrations::compute(const Molecule& molecule, const Hessian& hessian) {
    frequencies.clear();
    eigenvectors.clear();

    int n_atoms = molecule.get_num_atoms();
    int dim = molecule.get_num_coordinates();

    if (hessian.get_size() != dim) {
        throw std::runtime_error(
            "Vibrations::compute: Hessian size does not match 3N coordinates."
        );
    }

    std::vector<double> masses = molecule.get_masses();  // amu

    // Load Hessian into Eigen matrix:
    
    Eigen::MatrixXd H(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            H(i, j) = hessian.get_value(i, j);
        }
    }

    // Make sure Hessian is symmetric:
    
    H = 0.5 * (H + H.transpose());

    // Mass-weight Hessian using atomic-unit masses:
    
    Eigen::MatrixXd MW(dim, dim);
    for (int i = 0; i < dim; i++) {
        double mi = masses[i / 3] * AMU_TO_AU;

        for (int j = 0; j < dim; j++) {
            double mj = masses[j / 3] * AMU_TO_AU;
            MW(i, j) = H(i, j) / std::sqrt(mi * mj);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(MW);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Vibrations::compute: diagonalization failed.");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvecs = solver.eigenvectors();

    struct Mode {
        double freq_cm;
        Eigen::VectorXd vec;
    };

    std::vector<Mode> all_modes;

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);

        double freq = 0.0;
        if (lambda >= 0.0) {
            freq = std::sqrt(lambda) * AU_TO_CM;
        } else {
            freq = -std::sqrt(-lambda) * AU_TO_CM;
        }

        all_modes.push_back({freq, eigenvecs.col(i)});
    }

    std::sort(all_modes.begin(), all_modes.end(),
              [](const Mode& a, const Mode& b) {
                  return a.freq_cm < b.freq_cm;
              });

    int expected_modes = is_linear_molecule(molecule)
                         ? 3 * n_atoms - 5
                         : 3 * n_atoms - 6;

    std::vector<Mode> real_modes;
    for (const auto& mode : all_modes) {
        if (std::abs(mode.freq_cm) > FREQ_ZERO_THRESHOLD) {
            real_modes.push_back(mode);
        }
    }

    // If there are extra modes, remove the lowest-frequency ones:
    
    if (static_cast<int>(real_modes.size()) > expected_modes) {
        int extra = static_cast<int>(real_modes.size()) - expected_modes;
        real_modes.erase(real_modes.begin(), real_modes.begin() + extra);
    }

    for (const auto& mode : real_modes) {
        frequencies.push_back(mode.freq_cm);
        eigenvectors.push_back(mode.vec);
    }

    std::cout << "Found " << frequencies.size()
              << " vibrational modes (expected " << expected_modes << ").\n";
}

void Vibrations::print_frequencies() const {
    std::cout << "\n Vibrational Frequencies \n";

    for (size_t i = 0; i < frequencies.size(); i++) {
        std::cout << "  Mode " << i + 1
                  << ": " << frequencies[i] << " cm^-1\n";
    }

    std::cout << "----------------------------------------------\n";
    std::cout << "Reference values (Psi4/experimental value):\n";
    std::cout << "  H2:  ~4400 cm^-1  (stretch)\n";
    std::cout << "  HCl: ~2990 cm^-1  (stretch)\n";
    std::cout << "  H2O: ~1595, ~3657, ~3756 cm^-1\n";
}

void Vibrations::write_frequencies(const std::string& filename) const {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error(
            "Vibrations::write_frequencies: cannot write to " + filename
        );
    }

    file << "# Vibrational frequencies (cm^-1)\n";
    for (double f : frequencies) {
        file << f << "\n";
    }
}

std::vector<double> Vibrations::get_frequencies() const {
    return frequencies;
}

void Vibrations::write_normal_mode_xyz(
    const Molecule& molecule,
    const std::string& filename,
    int n_frames,
    double amplitude
) const {
    if (eigenvectors.empty()) {
        throw std::runtime_error("write_normal_mode_xyz: call compute() first.");
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error(
            "write_normal_mode_xyz: cannot write to " + filename
        );
    }

    int n_atoms = molecule.get_num_atoms();
    auto symbols = molecule.get_symbols();
    auto coords = molecule.get_coordinates();

    for (size_t m = 0; m < eigenvectors.size(); m++) {
        const Eigen::VectorXd& evec = eigenvectors[m];

        double norm = evec.norm();
        if (norm < 1.0e-12) {
            continue;
        }

        for (int frame = 0; frame < n_frames; frame++) {
            double phase = static_cast<double>(frame) /
                           static_cast<double>(n_frames);
            double disp = amplitude * std::sin(2.0 * M_PI * phase);

            file << n_atoms << "\n";
            file << "Mode " << m + 1
                 << "  nu=" << frequencies[m] << "_cm-1"
                 << "  frame=" << frame << "\n";

            for (int A = 0; A < n_atoms; A++) {
                double x = coords[3 * A]     + disp * evec(3 * A)     / norm;
                double y = coords[3 * A + 1] + disp * evec(3 * A + 1) / norm;
                double z = coords[3 * A + 2] + disp * evec(3 * A + 2) / norm;

                file << symbols[A] << "  "
                     << x << "  " << y << "  " << z << "\n";
            }
        }
    }

    std::cout << "Normal mode animation written to: " << filename << "\n";
}
