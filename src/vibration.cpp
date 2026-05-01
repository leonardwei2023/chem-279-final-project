#include "vibration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

// -----------------------------------------------------------------------
// Unit conversion: Hessian (Hartree/Bohr^2 / amu) -> frequency (cm^-1)
//
// From the proposal and Lecture 15:
//   omega_i = sqrt(lambda_i)        [angular frequency in atomic units]
//   nu_i    = omega_i / (2*pi*c)    [convert to wavenumbers, cm^-1]
//
// Constants:
//   1 Hartree = 4.3597447e-18 J
//   1 Bohr    = 5.2917721e-11 m
//   1 amu     = 1.6605391e-27 kg
//   c         = 2.9979246e10  cm/s
//
// The combined conversion factor is:
//   AU_TO_CM = sqrt(Ha/Bohr^2/amu) / (2*pi*c) = 5140.48  cm^-1
//
// Note: finite_difference.cpp divides the raw Hessian by ang_to_bohr^2
// to convert from Hartree/Ang^2 to Hartree/Bohr^2.  After mass-weighting
// by amu, each eigenvalue lambda is in Hartree/(Bohr^2*amu), so
// AU_TO_CM = 5140.48 is the correct constant.
// -----------------------------------------------------------------------
static constexpr double AU_TO_CM = 5140.48;

// Modes below this threshold (cm^-1) are treated as zero modes
// (translations + rotations) and discarded.
static constexpr double FREQ_ZERO_THRESHOLD = 10.0;

bool Vibrations::is_linear_molecule(const Molecule& molecule) const {
    return molecule.get_num_atoms() == 2;
}

void Vibrations::compute(const Molecule& molecule, const Hessian& hessian) {
    frequencies.clear();
    eigenvectors.clear();

    int n_atoms = molecule.get_num_atoms();
    int dim     = molecule.get_num_coordinates();  // 3 * n_atoms

    if (hessian.get_size() != dim) {
        throw std::runtime_error(
            "Vibrations::compute: Hessian size (" +
            std::to_string(hessian.get_size()) +
            ") != 3*N_atoms (" +
            std::to_string(dim) + ")."
        );
    }

    std::vector<double> masses = molecule.get_masses();  // amu

    // ------------------------------------------------------------------
    // Step 1: Load Hessian into Eigen and symmetrize
    // Finite differences can produce small asymmetries; symmetrizing
    // ensures SelfAdjointEigenSolver receives an exactly symmetric matrix.
    // ------------------------------------------------------------------
    Eigen::MatrixXd H(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            H(i, j) = hessian.get_value(i, j);
        }
    }
    H = 0.5 * (H + H.transpose());

    // ------------------------------------------------------------------
    // Step 2: Mass-weight the Hessian
    //   H'_{ij} = H_{ij} / sqrt(m_i * m_j)
    // (Proposal section b-iii; Lecture 15 vibrations slide)
    // ------------------------------------------------------------------
    Eigen::MatrixXd MW(dim, dim);
    for (int i = 0; i < dim; i++) {
        double mi = masses[i / 3];
        for (int j = 0; j < dim; j++) {
            double mj = masses[j / 3];
            MW(i, j) = H(i, j) / std::sqrt(mi * mj);
        }
    }

    // ------------------------------------------------------------------
    // Step 3: Diagonalize H' q = lambda q
    // SelfAdjointEigenSolver returns eigenvalues in ascending order.
    // ------------------------------------------------------------------
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(MW);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Vibrations::compute: Eigen diagonalization failed."
        );
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvecs   = solver.eigenvectors();

    // ------------------------------------------------------------------
    // Step 4: Convert eigenvalues to frequencies (cm^-1)
    //   nu_i = sqrt(lambda_i) * AU_TO_CM   for lambda > 0 (real mode)
    //   nu_i = -sqrt(|lambda_i|) * AU_TO_CM for lambda < 0 (imaginary)
    // (Proposal section b-iv)
    // ------------------------------------------------------------------
    struct Mode {
        double freq_cm;
        Eigen::VectorXd vec;
    };
    std::vector<Mode> all_modes;

    for (int i = 0; i < eigenvalues.size(); i++) {
        double lambda = eigenvalues(i);
        double freq   = (lambda >= 0.0)
                        ?  std::sqrt(lambda)  * AU_TO_CM
                        : -std::sqrt(-lambda) * AU_TO_CM;
        all_modes.push_back({freq, eigenvecs.col(i)});
    }

    // Sort ascending (negatives first, zero modes, then real vibrations)
    std::sort(all_modes.begin(), all_modes.end(),
              [](const Mode& a, const Mode& b) {
                  return a.freq_cm < b.freq_cm;
              });

    // ------------------------------------------------------------------
    // Step 5: Remove near-zero modes (translations + rotations)
    //   Linear molecule:    5 zero modes  -> 3N-5 vibrational modes
    //   Nonlinear molecule: 6 zero modes  -> 3N-6 vibrational modes
    // ------------------------------------------------------------------
    int expected_modes = is_linear_molecule(molecule)
                         ? 3 * n_atoms - 5
                         : 3 * n_atoms - 6;

    std::vector<Mode> real_modes;
    for (auto& m : all_modes) {
        if (std::abs(m.freq_cm) > FREQ_ZERO_THRESHOLD) {
            real_modes.push_back(m);
        }
    }

    // If still more than expected, drop lowest-|freq| extras
    if (static_cast<int>(real_modes.size()) > expected_modes) {
        real_modes.erase(
            real_modes.begin(),
            real_modes.begin() +
                (static_cast<int>(real_modes.size()) - expected_modes)
        );
    }

    for (auto& m : real_modes) {
        frequencies.push_back(m.freq_cm);
        eigenvectors.push_back(m.vec);
    }

    std::cout << "Found " << frequencies.size()
              << " vibrational modes (expected " << expected_modes << ").\n";
}

void Vibrations::print_frequencies() const {
    std::cout << "\n--- Vibrational Frequencies ---\n";
    for (size_t i = 0; i < frequencies.size(); i++) {
        std::cout << "  Mode " << i + 1
                  << ": " << frequencies[i] << " cm^-1\n";
    }
    std::cout << "-------------------------------\n";
    std::cout << "Reference values (Psi4/experimental):\n";
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

// -----------------------------------------------------------------------
// Write normal mode animation as a multi-frame .xyz file.
// Atoms are displaced along the normal mode eigenvector as a sine wave.
// Open the output in Avogadro or VMD for visualization.
// -----------------------------------------------------------------------
void Vibrations::write_normal_mode_xyz(
    const Molecule& molecule,
    const std::string& filename,
    int    n_frames,
    double amplitude
) const {
    if (eigenvectors.empty()) {
        throw std::runtime_error(
            "write_normal_mode_xyz: call compute() first."
        );
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error(
            "write_normal_mode_xyz: cannot write to " + filename
        );
    }

    int  n_atoms  = molecule.get_num_atoms();
    auto symbols  = molecule.get_symbols();
    auto coords   = molecule.get_coordinates();  // Angstrom, flat

    for (size_t m = 0; m < eigenvectors.size(); m++) {
        const Eigen::VectorXd& evec = eigenvectors[m];
        double norm = evec.norm();
        if (norm < 1.0e-12) continue;

        for (int frame = 0; frame < n_frames; frame++) {
            double phase = static_cast<double>(frame) /
                           static_cast<double>(n_frames);
            double disp  = amplitude * std::sin(2.0 * M_PI * phase);

            file << n_atoms << "\n";
            file << "Mode " << m + 1
                 << "  nu=" << frequencies[m] << "_cm-1"
                 << "  frame=" << frame << "\n";

            for (int A = 0; A < n_atoms; A++) {
                double x = coords[3*A]   + disp * evec(3*A)   / norm;
                double y = coords[3*A+1] + disp * evec(3*A+1) / norm;
                double z = coords[3*A+2] + disp * evec(3*A+2) / norm;
                file << symbols[A] << "  "
                     << x << "  " << y << "  " << z << "\n";
            }
        }
    }

    std::cout << "Normal mode animation written to: " << filename << "\n";
    std::cout << "  " << eigenvectors.size()
              << " modes, " << n_frames << " frames each.\n";
    std::cout << "  Open with Avogadro or VMD.\n";
}
