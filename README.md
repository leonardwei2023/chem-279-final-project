# CNDO/2 Vibrational Frequencies and Dipole Moments

**CHEM 279 — Numerical Algorithms Applied to Computational Quantum Chemistry**  
**Final Project — Spring 2026**

David Houshangi · Leonard Ming Wei  
University of California, Berkeley

---

## Overview

This project extends a CNDO/2 semi-empirical quantum chemistry implementation to compute two molecular properties beyond total energy:

1. **Vibrational frequencies** via a finite-difference Hessian
2. **Molecular dipole moments** from the CNDO/2 density matrix

The theory follows Lectures 13–17 of CHEM 279.

All equations are cited inline in the source code.

---

## Theory

### CNDO/2 SCF (Lecture 13)

The Complete Neglect of Differential Overlap (CNDO/2) model simplifies ab-initio Hartree-Fock theory by applying the Zero Differential Overlap (ZDO) approximation. Only valence electrons are treated explicitly:

- **H**: 1s orbital (1 basis function)
- **C, N, O, F, Cl**: 2s, 2px, 2py, 2pz (4 basis functions, no core 1s)

The Fock matrix elements are (Lec 13, slide 7):

```
Diagonal:
  F_μμ = -½(I_μ + A_μ) + [P_AA^tot - Z_A - (P_μμ^α - ½)]·γ_AA + Σ_{B≠A} [P_BB^tot - Z_B]·γ_AB

Off-diagonal (μ on A, ν on B, A≠B):
  F_μν = -½(β_A + β_B)·S_μν - P_μν^α · γ_AB
```

where `P^α = P^total / 2` for closed-shell molecules.

### Gamma Integrals — Mataga-Nishimoto (Lecture 14)

```
γ_AA = I_s - A_s                              (on-site)
γ_AB = 1 / (R_AB + 2 / (γ_AA + γ_BB))        (two-center)
```

### Overlap Integrals (Lectures 7–8, 10)

Computed analytically using STO-3G contracted Gaussian basis functions (Hehre, Stewart & Pople 1969).

### Finite-Difference Hessian (Lectures 2, 15)

Second derivatives of total energy using central differences:

```
Diagonal:    H_ii = [E(x+h) - 2E(x) + E(x-h)] / h²
Off-diagonal: H_ij = [E(++)-E(+-)-E(-+)+E(--)] / (4h²)
```

Raw values in Hartree/Å² are converted to Hartree/Bohr² by dividing by `(1.8897)²`.

### Vibrational Frequencies (Lecture 15, Proposal)

1. Mass-weight the Hessian: `H'_ij = H_ij / √(m_i · m_j)`
2. Diagonalize: `H' q = λ q`
3. Convert eigenvalues to frequencies: `ν_i = √λ_i × 5140.49 cm⁻¹`

The conversion factor `5140.49 cm⁻¹` comes from:

```
ν = √(Ha/Bohr²/amu) / (2π·c)
  = √(4.3597e-18 J / (2.8003e-21 m²) / (1.6605e-27 kg)) / (2π × 2.9979e10 cm/s)
  = 5140.49 cm⁻¹
```

Linear molecules have `3N−5` vibrational modes; nonlinear have `3N−6`.

### Dipole Moment (Lecture 16)

Under CNDO/2 ZDO, the SCF dipole moment (Lec 16, slide 10) reduces to:

```
μ_j = Σ_A (Z_A - P_AA) · R_A^j    [in e·Å]
|μ|  = √(μ_x² + μ_y² + μ_z²)     × 4.80320  [Debye]
```

where `P_AA = Σ_{μ on A} P_μμ` is the electron population on atom A.

---

## Atomic Parameters

Standard CNDO/2 parameters from Beveridge & Pople:

| Element | Z_val | I_s (eV) | A_s (eV) | I_p (eV) | A_p (eV) | β (eV) | ζ_s   | ζ_p   |
|---------|-------|----------|----------|----------|----------|--------|-------|-------|
| H       | 1     | 13.06    | 3.69     | —        | —        | −9.0   | 1.200 | —     |
| C       | 4     | 14.051   | 5.572    | 5.572    | 2.275    | −21.0  | 1.625 | 1.625 |
| N       | 5     | 19.316   | 7.275    | 7.275    | 2.862    | −25.0  | 1.950 | 1.950 |
| O       | 6     | 25.390   | 9.111    | 9.111    | 3.524    | −31.0  | 2.275 | 2.275 |
| F       | 7     | 32.272   | 11.080   | 11.080   | 4.270    | −39.0  | 2.425 | 2.425 |
| Cl      | 7     | 24.350   | 3.890    | 3.890    | 1.430    | −15.0  | 3.500 | 2.033 |

---

## Repository Structure

```
cndo2_vibrational/
├── CMakeLists.txt
├── README.md
├── .gitignore
│
├── src/
│   ├── main.cpp               # Program entry point, mode dispatch
│   ├── molecule.cpp           # Reads/writes .xyz files, stores coordinates and masses
│   ├── cndo_energy.cpp        # Full CNDO/2 SCF energy (called by CNDOEngine)
│   ├── cndo_engine.cpp        # Runs cndo_energy binary via shell, extracts energy
│   ├── finite_difference.cpp  # Builds Hessian via central differences
│   ├── hessian.cpp            # Hessian matrix: storage, read/write
│   ├── vibration.cpp          # Mass-weighted diagonalization → frequencies + modes
│   ├── dipole.cpp             # Dipole moment from CNDO/2 ZDO formula
│   └── validation.cpp         # Compares computed vs reference frequencies
│
├── include/
│   ├── molecule.h
│   ├── cndo_engine.h
│   ├── finite_difference.h
│   ├── hessian.h
│   ├── vibration.h
│   ├── dipole.h
│   └── validation.h
│
├── examples/
│   ├── h2.xyz                 # Hydrogen molecule (equilibrium geometry)
│   ├── hcl.xyz                # Hydrogen chloride
│   ├── h2o.xyz                # Water
│   ├── h2_hessian.dat         # Pre-computed Hessian for H2
│   ├── hcl_hessian.dat        # Pre-computed Hessian for HCl
│   └── h2o_hessian.dat        # Pre-computed Hessian for H2O
│
└── scripts/
    └── psi4_reference.py      # Generate reference frequencies with Psi4/HF/STO-3G
```

---

## Dependencies

- **C++17** or later
- **Eigen3** — header-only linear algebra library
  - Ubuntu/Debian: `sudo apt install libeigen3-dev`
  - macOS: `brew install eigen`
- **CMake** ≥ 3.10

---

## Build

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

This produces two executables in `build/`:
- `vibrational_frequency` — main program
- `cndo_energy` — CNDO/2 energy oracle (called internally by `vibrational_frequency`)

---

## Run

### 1. Vibrational frequencies from a pre-computed Hessian

```bash
./vibrational_frequency vibration ../examples/hcl.xyz ../examples/hcl_hessian.dat
```

Add `--animate` to write `normal_modes.xyz` for Avogadro/VMD:

```bash
./vibrational_frequency vibration ../examples/h2o.xyz ../examples/h2o_hessian.dat --animate
```

### 2. Compute Hessian via finite differences (requires CNDO/2 energy)

```bash
export CNDO_ENERGY_CMD="./cndo_energy"
./vibrational_frequency finite-diff ../examples/hcl.xyz hcl_computed.dat 0.01
```

### 3. Compute Hessian AND frequencies in one step

```bash
export CNDO_ENERGY_CMD="./cndo_energy"
./vibrational_frequency finite-diff-vib ../examples/h2.xyz h2_computed.dat 0.005 --animate
```

### 4. Dipole moment (neutral-atom approximation, placeholder)

```bash
./vibrational_frequency dipole ../examples/hcl.xyz
```

### 5. Dipole moment with SCF density matrix

Once you have per-atom electron populations from CNDO/2 (one value per atom, `sum_{μ on A} P_{μμ}`):

```bash
./vibrational_frequency dipole-scf ../examples/hcl.xyz p_diagonal.dat
```

### 6. Validate computed frequencies against reference

```bash
./vibrational_frequency validate computed_frequencies.dat ../examples/hcl_hessian.dat
```

---

## Expected Result Values

| Molecule | Mode        | Computed (cm⁻¹) | Experiment (cm⁻¹) |
|----------|-------------|-----------------|-------------------|
| H₂       | stretch     | ~4100–4400      | 4161              |
| HCl      | stretch     | ~2700–3100      | 2991              |
| H₂O      | bend        | ~1300–1700      | 1595              |
| H₂O      | sym stretch | ~3200–3700      | 3657              |
| H₂O      | asym stretch| ~3300–3800      | 3756              |

CNDO/2 is a semi-empirical method; errors of 10–20% relative to experiment are expected and consistent with the method's known limitations (Lec 13, Lec 17).

---

## Usage Notes

- **Step size**: `h = 0.005` Å works well for CNDO/2 (balances truncation and roundoff error per Lec 2)
- **Temporary files**: `cndo_energy` writes and deletes `tmp_cndo_displaced_N.xyz` during finite differences
- **Open-shell molecules**: not supported (closed-shell RHF only)
- **Supported elements**: H, C, N, O, F, Cl

---

## File Format: p_diagonal.dat

For the `dipole-scf` mode, provide a plain text file with one number per line — the total electron population on each atom, in the same order as the `.xyz` file:

```
# p_diagonal.dat for HCl: H Cl
0.72
7.28
```

This is `sum_{μ on A} P_{μμ}` from the converged CNDO/2 density matrix.

---

## Generating Reference Frequencies with Psi4

```bash
python3 scripts/psi4_reference.py examples/hcl.xyz hcl_reference.dat
./vibrational_frequency validate computed_frequencies.dat hcl_reference.dat
```

---

## Authors

David Houshangi — davidhoushangi@berkeley.edu  
Leonard Ming Wei — dranoelmi@berkeley.edu

CHEM 279, Spring 2026  
Instructor: Dr. Mayank Agrawal

---

## Contributions

**David Houshangi**: CNDO/2 SCF implementation (`cndo_energy.cpp`), overlap integrals, gamma integrals, Fock matrix construction, SCF loop.

**Leonard Ming Wei**: Finite-difference Hessian (`finite_difference.cpp`), vibrational frequency analysis (`vibration.cpp`), dipole moment (`dipole.cpp`), validation framework, normal mode animation output.
