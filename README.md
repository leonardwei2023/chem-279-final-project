# CNDO/2 Vibrational Frequencies and Dipole Moments

**CHEM 279 — Numerical Algorithms Applied to Computational Quantum Chemistry**
**Final Project — Spring 2026**

David Houshangi · Leonard Ming Wei
University of California, Berkeley

---

## Overview

This project implements a **CNDO/2 semi-empirical quantum chemistry model** to compute:

1. **Vibrational frequencies** via a finite-difference Hessian
2. **Molecular dipole moments** from the **SCF density matrix**

All computations follow Lectures 13–17 of CHEM 279.

---

## Theory

### CNDO/2 SCF (Lecture 13)

CNDO/2 applies the **Zero Differential Overlap (ZDO)** approximation and treats only valence electrons.

The Fock matrix:

Diagonal:
Fμμ = −½(Iμ + Aμ) + [PAA − ZA − (Pμμ^α − ½)] γAA + ΣB≠A (PBB − ZB) γAB

Off-diagonal:
Fμν = −½(βA + βB) Sμν − Pμν^α γAB

---

### Density Matrix

Closed-shell density:

P = 2 C_occ C_occᵀ

Atomic populations:

PAA = Σ (μ on atom A) Pμμ

---

### Dipole Moment (Lecture 16)

Dipole is computed from SCF density:

μ = ΣA (ZA − PAA) RA   (e·Å)

Converted to Debye:

|μ| = √(μx² + μy² + μz²) × 4.80320

---

### Finite-Difference Hessian

Hessian elements:

Hii = [E(x+h) − 2E(x) + E(x−h)] / h²
Hij = [E(++) − E(+-) − E(-+) + E(--)] / (4h²)

Converted from Å → Bohr units.

---

### Vibrational Frequencies

Mass-weighted Hessian:

H'ij = Hij / √(mi mj)

Eigenvalues:

ν = √λ × 5140.49 cm⁻¹

---

## Build

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

---

## Run

### Vibrations (from Hessian)

```bash
./vibrational_frequency vibration ../examples/hcl.xyz ../examples/hcl_hessian.dat
```

---

### Finite Difference + Frequencies

```bash
export CNDO_ENERGY_CMD="./cndo_energy"
./vibrational_frequency finite-diff-vib ../examples/h2.xyz h2.dat 0.005
```

---

### Dipole Moment (FULL SCF — CORRECT)

```bash
./vibrational_frequency dipole ../examples/h2o.xyz
```

This automatically:

1. Runs CNDO/2 SCF
2. Generates `p_diagonal.dat`
3. Computes dipole using SCF density

---

## Example Output

```text
SCF energy = -12.2475 Hartree

Dipole Moment (CNDO/2 SCF Density)
|mu| = 0.45 Debye
```

---

## Expected Values

| Molecule | Property | CNDO/2          | Experimental |
| -------- | -------- | --------------- | ------------ |
| H2       | stretch  | ~4100–4400 cm⁻¹ | 4161         |
| HCl      | stretch  | ~2700–3100 cm⁻¹ | 2991         |
| H2O      | dipole   | ~0.4–1.0 Debye  | 1.85         |

---

## Notes

* CNDO/2 is a **semi-empirical model**, so deviations from experiment are expected
* No artificial charge fitting is used
* Dipole is computed **directly from SCF density matrix**
* Fully consistent with lecture formulation

---

## Authors

David Houshangi
Leonard Ming Wei

CHEM 279 — Spring 2026
