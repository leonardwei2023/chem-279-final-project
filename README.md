# CNDO/2 Vibrational Frequencies and Dipole Moments

**CHEM 279 — Final Project**

David Houshangi · Leonard Ming Wei

UC Berkeley

---

## Overview

This project implements a small quantum chemistry workflow using the **CNDO/2 semi-empirical method**.

The program can:

* Compute **vibrational frequencies** from a Hessian matrix
* Generate a **finite-difference Hessian** using CNDO/2 energies
* Compute **dipole moments** using the SCF density matrix

Everything is written in C++ and built with CMake.

---

## What This Project Does

### 1. CNDO/2 SCF Energy

The code performs a self-consistent field (SCF) calculation to get:

* Molecular energy
* Density matrix
* Atomic electron populations (written to `p_diagonal.dat`)

These populations are then used for dipole calculations.

---

### 2. Vibrational Frequencies

Frequencies are computed by:

1. Building a **Hessian matrix**
2. Converting to a **mass-weighted Hessian**
3. Solving for eigenvalues

The output is in:

```text
cm⁻¹
```

---

### 3. Finite-Difference Hessian

Instead of reading a Hessian from file, the program can compute it numerically:

* Displace coordinates by a small step
* Recompute energy using CNDO/2
* Build second derivatives

This is slower (~20–30 sec per molecule) but more realistic.

---

### 4. Dipole Moment (SCF-Based)

The dipole moment is computed from the electron density:

```text
μ = ΣA (Z_A − P_AA) R_A
```

Where:

* `Z_A` = nuclear charge
* `P_AA` = electron population on atom A
* `R_A` = atomic position

The result is converted to:

```text
Debye
```

---

## Project Structure

```bash
.
├── include/        # Header files
├── src/            # Source files
├── input/          # Test molecules (.xyz, .dat)
├── examples/       # Example runs
├── build/          # Build directory
├── CMakeLists.txt
└── README.md
```

---

## Build Instructions

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

---

## How to Run

### Vibrational Frequencies (from Hessian file)

```bash
./vibrational_frequency vibration ../input/hcl.xyz ../input/hcl_hessian.dat
```

---

### Finite Difference + Frequencies

```bash
export CNDO_ENERGY_CMD="./cndo_energy"

./vibrational_frequency finite-diff-vib \
    ../input/h2.xyz h2_hessian.dat 0.005
```

---

### Dipole Moment (Recommended)

```bash
./vibrational_frequency dipole ../input/h2o.xyz
```

This will:

1. Run CNDO/2 SCF
2. Generate `p_diagonal.dat`
3. Compute the dipole moment

---

### Dipole Using Existing SCF File

```bash
./vibrational_frequency dipole-scf ../input/h2o.xyz p_diagonal.dat
```

---

## Example Output

```text
SCF energy = -12.2475 Hartree

Dipole Moment
mu_x = 0
mu_y = 0
mu_z = -0.450 Debye
|mu| = 0.450 Debye
```

---

## Expected Results

| Molecule | Property          | Typical Output  |
| -------- | ----------------- | --------------- |
| H2       | Stretch frequency | ~4100–4400 cm⁻¹ |
| HCl      | Stretch frequency | ~2700–3100 cm⁻¹ |
| H2O      | Dipole moment     | ~0.4–1.0 Debye  |

Experimental values are usually higher because this is a simplified model.

---

## Notes

* Only **closed-shell systems** are supported
* Uses **valence electrons only**
* Finite-difference step size affects accuracy
* No geometry optimization is performed
* Dipole moment comes directly from SCF density (no manual tuning)

---

## Output Files

| File                       | Description                 |
| -------------------------- | --------------------------- |
| `p_diagonal.dat`           | Atomic electron populations |
| `computed_frequencies.dat` | Vibrational frequencies     |
| `normal_modes.xyz`         | Animation file (optional)   |
| `dipole_moment.dat`        | Dipole vector               |

---

## Tips

* Use a smaller step size (e.g., `0.005`) for better Hessians
* Always run dipole through SCF (not manual input)
* Check that atom ordering matches between files

---

## Authors

David Houshangi

Leonard Ming Wei

Spring 2026
