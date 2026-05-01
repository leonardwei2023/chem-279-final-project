#include "finite_difference.h"
#include <iostream>

FiniteDifference::FiniteDifference(double step) {
    step_size = step;
}

Hessian FiniteDifference::compute_hessian(const Molecule& molecule, CNDOEngine& engine) {
    int size = molecule.get_num_coordinates();

    Hessian hessian;
    hessian.resize(size);

    const double ang_to_bohr = 1.889726124565062;
    double step_bohr = step_size * ang_to_bohr;

    std::cout << "Computing finite-difference Hessian\n";
    std::cout << "Coordinates: " << size << "\n";
    std::cout << "Step size: " << step_size << " Angstrom\n";
    std::cout << "Step size: " << step_bohr << " Bohr\n";

    double e0 = engine.compute_energy(molecule);

    for (int i = 0; i < size; i++) {
        for (int j = i; j < size; j++) {
            double value = 0.0;

            if (i == j) {
                Molecule mol_plus = molecule;
                Molecule mol_minus = molecule;

                mol_plus.set_coordinate(i, mol_plus.get_coordinate(i) + step_size);
                mol_minus.set_coordinate(i, mol_minus.get_coordinate(i) - step_size);

                double e_plus = engine.compute_energy(mol_plus);
                double e_minus = engine.compute_energy(mol_minus);

                value = (e_plus - 2.0 * e0 + e_minus) /
                        (step_bohr * step_bohr);
            } else {
                Molecule mol_pp = molecule;
                Molecule mol_pm = molecule;
                Molecule mol_mp = molecule;
                Molecule mol_mm = molecule;

                mol_pp.set_coordinate(i, mol_pp.get_coordinate(i) + step_size);
                mol_pp.set_coordinate(j, mol_pp.get_coordinate(j) + step_size);

                mol_pm.set_coordinate(i, mol_pm.get_coordinate(i) + step_size);
                mol_pm.set_coordinate(j, mol_pm.get_coordinate(j) - step_size);

                mol_mp.set_coordinate(i, mol_mp.get_coordinate(i) - step_size);
                mol_mp.set_coordinate(j, mol_mp.get_coordinate(j) + step_size);

                mol_mm.set_coordinate(i, mol_mm.get_coordinate(i) - step_size);
                mol_mm.set_coordinate(j, mol_mm.get_coordinate(j) - step_size);

                double e_pp = engine.compute_energy(mol_pp);
                double e_pm = engine.compute_energy(mol_pm);
                double e_mp = engine.compute_energy(mol_mp);
                double e_mm = engine.compute_energy(mol_mm);

                value = (e_pp - e_pm - e_mp + e_mm) /
                        (4.0 * step_bohr * step_bohr);
            }

            hessian.set_value(i, j, value);
            hessian.set_value(j, i, value);
        }

        std::cout << "Finished Hessian row " << i + 1
                  << " of " << size << "\n";
    }

    return hessian;
}
