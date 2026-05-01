#ifndef VIBRATION_H
#define VIBRATION_H

#include "hessian.h"
#include "molecule.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

/**
 * Vibrations: computes vibrational frequencies from the Hessian.
 *
 * Algorithm (Lecture 15):
 
 *   1. Symmetrize the Hessian
 *   2. Mass-weight:  H'_{ij} = H_{ij} / sqrt(m_i * m_j)
 *   3. Diagonalize:  H' q = lambda q
 *   4. Convert:      nu_i = sqrt(lambda_i) * AU_TO_CM  (cm^-1)
 *   5. Filter out 5 (linear) or 6 (nonlinear) zero modes
 */

class Vibrations {
public:
    // Compute frequencies and store normal mode eigenvectors

    void compute(const Molecule& molecule, const Hessian& hessian);

    // Print frequencies to stdout with reference values

    void print_frequencies() const;

    // Write frequencies (one per line) to file:

    void write_frequencies(const std::string& filename) const;

    /**
     * Write multi-frame .xyz animation for all normal modes.
     * Load in Avogadro or VMD to visualize atomic motion.
     *
     * @param molecule   The molecule
     * @param filename   Output .xyz file path
     * @param n_frames   Frames per mode (default 20)
     * @param amplitude  Max displacement in Angstrom (default 0.3)
     */

    void write_normal_mode_xyz(
        const Molecule& molecule,
        const std::string& filename,
        int    n_frames  = 20,
        double amplitude = 0.3
    ) const;

    std::vector<double> get_frequencies() const;

private:
    std::vector<double>          frequencies;
    std::vector<Eigen::VectorXd> eigenvectors;

    bool is_linear_molecule(const Molecule& molecule) const;
};

#endif // VIBRATION_H
