#include "vibrational_frequency.h"

#include <armadillo>
#include <cmath>
#include <iomanip>
#include <iostream>

// Temporary energy function.
// This is not real CNDO/2 yet. It is a simple harmonic energy model
// so we can test the vibrational-frequency pipeline first.
//
// E = 1/2 * k * (r - r0)^2


double compute_cndo2_total_energy(Molecule& mol)
{
    double energy = 0.0;
    int natoms = mol.num_atoms();

    for (int i = 0; i < natoms; i++) {
        for (int j = i + 1; j < natoms; j++) {
            arma::vec ri = mol.coordinates.row(i).t();
            arma::vec rj = mol.coordinates.row(j).t();

            double distance = arma::norm(ri - rj);

            double r0 = 1.0;
            double k = 0.5;

            energy += 0.5 * k * std::pow(distance - r0, 2);
        }
    }

    return energy;
}

// Unit conversion constants.
// Assumptions:
// energy   = Hartree
// distance = Bohr
// mass     = amu
// output   = cm^-1
const double HARTREE_TO_JOULE = 4.3597447222071e-18;
const double BOHR_TO_METER = 5.29177210903e-11;
const double AMU_TO_KG = 1.66053906660e-27;
const double SPEED_OF_LIGHT_CM = 2.99792458e10;

VibrationalFrequencyAnalyzer::VibrationalFrequencyAnalyzer(double step_size)
{
    h = step_size;
}

// Build Hessian using central finite differences.
//
// H_ij = d^2E / dxi dxj
//
// H_ij ≈ [E(++ ) - E(+-) - E(-+) + E(--)] / 4h^2
arma::mat VibrationalFrequencyAnalyzer::compute_hessian(Molecule mol)
{
    int natoms = mol.num_atoms();
    int ndim = 3 * natoms;

    arma::mat hessian(ndim, ndim, arma::fill::zeros);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Molecule mol_pp = mol;
            Molecule mol_pm = mol;
            Molecule mol_mp = mol;
            Molecule mol_mm = mol;

            int atom_i = i / 3;
            int coord_i = i % 3;

            int atom_j = j / 3;
            int coord_j = j % 3;

            mol_pp.coordinates(atom_i, coord_i) += h;
            mol_pp.coordinates(atom_j, coord_j) += h;

            mol_pm.coordinates(atom_i, coord_i) += h;
            mol_pm.coordinates(atom_j, coord_j) -= h;

            mol_mp.coordinates(atom_i, coord_i) -= h;
            mol_mp.coordinates(atom_j, coord_j) += h;

            mol_mm.coordinates(atom_i, coord_i) -= h;
            mol_mm.coordinates(atom_j, coord_j) -= h;

            double e_pp = compute_cndo2_total_energy(mol_pp);
            double e_pm = compute_cndo2_total_energy(mol_pm);
            double e_mp = compute_cndo2_total_energy(mol_mp);
            double e_mm = compute_cndo2_total_energy(mol_mm);

            hessian(i, j) =
                (e_pp - e_pm - e_mp + e_mm) / (4.0 * h * h);
        }
    }

    return hessian;
}

// Mass-weighted Hessian:
//
// H'_ij = H_ij / sqrt(m_i m_j)
arma::mat VibrationalFrequencyAnalyzer::mass_weight_hessian(
    const arma::mat& hessian,
    const Molecule& mol
)
{
    int natoms = mol.num_atoms();
    int ndim = 3 * natoms;

    arma::mat mass_weighted(ndim, ndim, arma::fill::zeros);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            int atom_i = i / 3;
            int atom_j = j / 3;

            double mass_i = mol.atomic_mass(atom_i);
            double mass_j = mol.atomic_mass(atom_j);

            mass_weighted(i, j) =
                hessian(i, j) / std::sqrt(mass_i * mass_j);
        }
    }

    return mass_weighted;
}

// Diagonalize mass-weighted Hessian:
//
// H' q = lambda q
//
// frequency = (1 / 2 pi c) * sqrt(lambda)
arma::vec VibrationalFrequencyAnalyzer::compute_frequencies(
    const arma::mat& mass_weighted_hessian
)
{
    arma::vec eigenvalues;
    arma::mat eigenvectors;

    arma::eig_sym(eigenvalues, eigenvectors, mass_weighted_hessian);

    arma::vec frequencies(eigenvalues.n_elem, arma::fill::zeros);

    double pi = std::acos(-1.0);

    double conversion =
        std::sqrt(
            HARTREE_TO_JOULE /
            (BOHR_TO_METER * BOHR_TO_METER * AMU_TO_KG)
        ) /
        (2.0 * pi * SPEED_OF_LIGHT_CM);

    for (arma::uword i = 0; i < eigenvalues.n_elem; i++) {
        double lambda = eigenvalues(i);

        if (lambda < 0.0) {
            frequencies(i) = -std::sqrt(std::abs(lambda)) * conversion;
        } else {
            frequencies(i) = std::sqrt(lambda) * conversion;
        }
    }

    return frequencies;
}

void VibrationalFrequencyAnalyzer::print_frequencies(
    const arma::vec& frequencies
)
{
    std::cout << "\nVibrational Frequency Analysis\n";
    std::cout << "------------------------------\n";
    std::cout << "Mode        Frequency (cm^-1)\n";

    for (arma::uword i = 0; i < frequencies.n_elem; i++) {
        std::cout << std::setw(4) << i + 1
                  << std::setw(20)
                  << std::fixed << std::setprecision(4)
                  << frequencies(i) << "\n";
    }

    std::cout << "\nNotes:\n";
    std::cout << "Linear molecules have 3N - 5 vibrational modes.\n";
    std::cout << "Nonlinear molecules have 3N - 6 vibrational modes.\n";
    std::cout << "This version uses a temporary harmonic energy model.\n";
}
