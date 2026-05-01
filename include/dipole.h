#ifndef DIPOLE_H
#define DIPOLE_H

/**
 * DipoleMoment: Computes the molecular electric dipole moment from CNDO/2.
 *
 * From Lecture 16 (Lec16, slide 10), the SCF dipole moment is:
 *
 *   mu = d_SCF = sum_{mu,nu} P_{mu,nu} * d_{mu,nu}
 *              = electronic part + nuclear part
 *
 * Where the dipole operator is:
 *   d_hat = -sum_i r_i  +  sum_A Z_A * R_A
 *            (electron)     (nuclear)
 *
 * In CNDO/2 with the ZDO approximation (S = I), off-diagonal dipole
 * integrals between orbitals on DIFFERENT atoms are zero.
 * Only on-site integrals survive:
 *   d_{mu,mu} = R_A  (position of atom A that owns orbital mu)
 *
 * So the electronic contribution reduces to:
 *   mu_elec = -sum_A [ sum_{mu on A} P_{mu,mu} ] * R_A
 *           = -sum_A q_A_elec * R_A
 *
 * And the full dipole is:
 *   mu = sum_A (Z_A - q_A_elec) * R_A
 *      = sum_A net_charge_A * R_A
 *
 * Units:
 *   - Coordinates in Angstrom, charges in electron units
 *   - 1 e * 1 Angstrom = 4.80320 Debye
 *   - Result reported in Debye
 */

#include "molecule.h"
#include <array>
#include <string>
#include <vector>

class DipoleMoment {
public:
    // Conversion factor: 1 e*Angstrom = 4.80320 Debye
    static constexpr double EA_TO_DEBYE = 4.80320;

    /**
     * Compute dipole moment using atomic partial charges.
     *
     * This uses the CNDO/2 ZDO-simplified formula:
     *   mu_j = sum_A (Z_A - P_AA) * R_A^j,  j in {x,y,z}
     *
     * @param molecule   The molecule (coordinates, symbols, valence Z)
     * @param p_diagonal Per-atom sum of diagonal density matrix elements.
     *                   p_diagonal[A] = sum_{mu on atom A} P_{mu,mu}
     *                   This is the number of electrons on atom A.
     *                   Pass an empty vector to use a nuclear-only estimate
     *                   (placeholder until SCF is connected).
     */
    void compute(
        const Molecule& molecule,
        const std::vector<double>& p_diagonal
    );

    // Print dipole components and magnitude to stdout
    void print() const;

    // Write dipole to file: x y z magnitude (Debye)
    void write(const std::string& filename) const;

    // Accessors
    std::array<double, 3> get_components_debye() const;
    double get_magnitude_debye() const;

private:
    // Dipole components in Debye
    std::array<double, 3> mu_debye = {0.0, 0.0, 0.0};
    double magnitude_debye = 0.0;

    // CNDO/2 valence electrons (Z_valence) for each element
    static int get_valence_electrons(const std::string& symbol);
};

#endif // DIPOLE_H
