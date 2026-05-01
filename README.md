# Calculation of Molecular Dipole Moments and Basic Vibrational Frequencies using CNDO/2

**CHEM 279 — Final Project**  

David Houshangi

Leonard Ming Wei  

UC Berkeley  

---

## Overview

This project implements a simple quantum chemistry workflow using the **CNDO/2 semi-empirical method**.

The program is able to:
- Compute **vibrational frequencies** using a finite-difference Hessian
- Generate a **Hessian matrix** directly from CNDO/2 energies
- Compute **dipole moments** using the SCF electron density

All calculations are performed internally (no precomputed data files are required).

---

## What This Project Does

### 1. CNDO/2 SCF Energy

The program performs a **self-consistent field (SCF)** calculation to obtain:

- Total molecular energy  
- Electron density (via atomic populations)  
- Atomic populations written to `p_diagonal.dat`  

These populations are later used for dipole moment calculations.

---

### 2. Vibrational Frequencies

Vibrational frequencies are computed through:

1. Building a **finite-difference Hessian**
2. Converting to a **mass-weighted Hessian**
3. Solving the eigenvalue problem

Final frequencies are reported in:

```
cm⁻¹
```

---

### 3. Finite-Difference Hessian

The Hessian is computed numerically:

- Each coordinate is displaced by a small step size  
- Energy is recomputed using CNDO/2  
- Second derivatives are assembled  

This creates a fully consistent pipeline:

```
Energy → Hessian → Vibrational Frequencies
```

---

### 4. Dipole Moment (SCF-Based)

The dipole moment is computed using:

```
μ = Σ_A (Z_A − P_AA) R_A
```

Where:

- `Z_A` = nuclear charge  
- `P_AA` = electron population on atom A  
- `R_A` = atomic position  

The final result is converted to:

```
Debye
```

---

## Project Structure

```
.
├── include/        # Header files
├── src/            # Source files
├── input/          # Molecule files (.xyz)
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

### 1. Vibrational Frequency (Recommended: H₂)

```bash
export CNDO_ENERGY_CMD="./cndo_energy"

./vibrational_frequency finite-diff ../input/h2.xyz h2_fd.dat 0.005
./vibrational_frequency vibration ../input/h2.xyz h2_fd.dat
```

---

### 2. Dipole Moment Calculations

#### Water (H₂O)
```bash
./vibrational_frequency dipole ../input/h2o.xyz
```

#### Hydrogen Chloride (HCl)
```bash
./vibrational_frequency dipole ../input/hcl.xyz
```

#### Ammonia (NH₃)
```bash
./vibrational_frequency dipole ../input/nh3.xyz
```

---

## Example Output

```
SCF energy = -12.2475 Hartree

Dipole Moment
mu_x = 0
mu_y = 0
mu_z = -0.450 Debye
|mu| = 0.450 Debye
```

---

## Results Summary

| Molecule | Property           | Result (This Work) | Experimental |
|----------|-------------------|--------------------|--------------|
| H₂       | Stretch frequency | ~4404 cm⁻¹         | ~4400 cm⁻¹   |
| H₂O      | Dipole moment     | ~0.45 Debye        | 1.85 Debye   |
| HCl      | Dipole moment     | ~0.24 Debye        | 1.08 Debye   |
| NH₃      | Dipole moment     | ~0.53 Debye        | 1.47 Debye   |

---

## Notes

- Only **closed-shell systems** are supported  
- Uses **valence electrons only (CNDO/2 approximation)**  
- No geometry optimization is performed  
- Finite-difference step size strongly affects accuracy  
- Dipole moments are computed directly from SCF populations  
- Vibrational analysis is most reliable for small systems (e.g., H₂)

---

## Output Files

| File               | Description                      |
|--------------------|----------------------------------|
| `p_diagonal.dat`   | Atomic electron populations      |
| `*_fd.dat`         | Generated Hessian matrices       |
| `normal_modes.xyz` | Optional vibrational animation   |

---

## Tips

- Use step size **0.003–0.005 Å** for stable Hessians  
- Always generate Hessians using finite differences  
- Ensure `.xyz` files are formatted correctly  
- Compare trends (not exact values) with experiment  

---

## Authors

David Houshangi  
Leonard Ming Wei  

Spring 2026
