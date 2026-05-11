# Qualitative Prediction of Vibrational Frequencies and Dipole Moments in Drug-Relevant Molecular Fragments Using CNDO/2 

> How Functional Groups Influence Polarity and Vibrational Signatures in Pharmaceutically Relevant Molecular Fragments.
- Authors: David Houshangi, Leonard Ming Wei
- CHEM 279 – Numerical Algorithms in Computational Quantum Chemistry 
- Methods: CNDO/2, SCF, Finite-Difference Hessian, Vibrational Analysis, Dipole Moments

## Project Overview

This project implements a modular computational quantum chemistry workflow using the CNDO/2 semiempirical method to study vibrational frequencies and dipole moments in small molecular systems relevant to common pharmaceutical functional groups.

The project includes:

- CNDO/2 self-consistent field (SCF) energy calculations
- Finite-difference Hessian generation
- Vibrational frequency calculations
- Dipole moment calculations from SCF electron populations
- Validation with Experimental Reference Data

Computed vibrational frequencies and dipole moments were compared against experimental reference values obtained from:

- NIST Chemistry WebBook
- NIST CCCBDB Database

The project focuses on reproducing qualitative molecular trends rather than exact quantitative agreement.
- Visualization of molecular vibrations and normal modes

The primary goal of this project is to evaluate whether CNDO/2 can reproduce qualitative trends in molecular polarity and vibrational behavior across chemically relevant functional groups.

Small molecular systems representing common functional groups were analyzed, including:

- H₂
- HCl
- H₂O
- NH₃
- Methanol
- Formaldehyde
- Acetaldehyde

The project focuses on understanding how functional groups influence molecular polarity, infrared-active vibrational modes, and spectroscopic behavior using an efficient semiempirical quantum chemistry framework.


> molecules for vibrational frequency calculations:
- H2
- HCl
- H2O

> molecules for dipole moment calculations:
- H2
- HCl
- H2O
- NH3
- Methanol
- Formaldehyde
- Acetaldehyde


# Theory

## 1. CNDO/2 Self-Consistent Field Method

The CNDO/2 method is a semiempirical quantum chemistry approximation that simplifies electronic integrals using the Complete Neglect of Differential Overlap (CNDO) approximation.

The SCF procedure iteratively computes:

- molecular electronic energy
- electron density
- atomic electron populations

Atomic populations are written to:

p_diagonal.dat

These populations are later used in dipole moment calculations.

---

## 2. Vibrational Frequencies

Vibrational frequencies are computed using a finite-difference Hessian approach.

The procedure consists of:

1. Displacing atomic coordinates by a small step size
2. Recomputing molecular energy using CNDO/2
3. Constructing the Hessian matrix from second derivatives
4. Building the mass-weighted Hessian
5. Solving the eigenvalue problem

Final frequencies are reported in: cm⁻¹

The implemented workflow creates a fully self-consistent pipeline:

Energy → Hessian → Vibrational Frequencies

---

## 3. Dipole Moment Calculation

The dipole moment is computed using:

μ = Σ_A (Z_A − P_AA) R_A

Where:

- Z_A = nuclear charge on atom A
- P_AA = SCF electron population on atom A
- R_A = atomic position vector

The final dipole moment is reported in: Debye

---

# Project Structure

```text
chem-279-final-project/
│
├── include/              # Header files
├── src/                  # Source files
├── input/                # Molecular geometry files (.xyz)
├── visualization/        # Visualization and plotting tools
├── output/               # Generated Hessians, frequencies, and outputs
├── report/               # Final report PDF and supporting files
├── slides/               # Presentation slides
│
├── CMakeLists.txt
├── Dockerfile
├── README.md
├── rubric.md
└── .gitignore
```

# Build Instructions

## Local Build

```bash
mkdir build
cd build

cmake ..
cmake --build .
```

# Running the Program

## Run Everything

From the main project directory:

```bash
chmod +x run.sh
./run.sh
```
## Get the out put for ./run.sh
```bash
chmod +x run.sh
./run.sh output.txt
```
All calculations are executed from the build directory.

Example:

```bash
cd build
./vibrational_frequency dipole ../input/h2o.xyz
# Vibrational Frequency Calculations
```

## Generate Hessian Matrices

### H₂
```
./vibrational_frequency finite-diff ../input/h2.xyz h2_fd.dat 0.005
```
## HCl
```
./vibrational_frequency finite-diff ../input/hcl.xyz hcl_fd.dat 0.005
```
## H2O
```
./vibrational_frequency finite-diff ../input/h2o.xyz h2o_fd.dat 0.005
```


# Compute Vibrational Frequencies from Hessians

## H₂
```
./vibrational_frequency vibration ../input/h2.xyz h2_fd.dat
```

## HCl
```
./vibrational_frequency vibration ../input/hcl.xyz hcl_fd.dat
```

## H2O
```
./vibrational_frequency vibration ../input/h2o.xyz h2o_fd.dat
```


# Compute Vibrational Frequencies Directly from Finite Differences

## H₂
```
./vibrational_frequency finite-diff-vib ../input/h2.xyz h2_fd.dat 0.005
```

## HCl
```
./vibrational_frequency finite-diff-vib ../input/hcl.xyz hcl_fd.dat 0.005
```

## H2O
```
./vibrational_frequency finite-diff-vib ../input/h2o.xyz h2o_fd.dat 0.005
```


# Dipole Moment Calculations

Water (H₂O)
```
./vibrational_frequency dipole ../input/h2o.xyz
```

Hydrogen Chloride (HCl)
```
./vibrational_frequency dipole ../input/hcl.xyz
```

Ammonia (NH₃)
```
./vibrational_frequency dipole ../input/nh3.xyz
```

Methanol
```
./vibrational_frequency dipole ../input/methanol.xyz
```

Formaldehyde
```
./vibrational_frequency dipole ../input/formaldehyde.xyz
```

Acetaldehyde
```
./vibrational_frequency dipole ../input/acetaldehyde.xyz
```

## Example Output

```
SCF energy = -12.2475 Hartree

Dipole Moment
mu_x = 0.000
mu_y = 0.000
mu_z = -0.450 Debye

|mu| = 0.450 Debye
```

# Input Molecules

The `input/` directory contains molecular geometry files in XYZ format used for vibrational frequency and dipole moment calculations.

| File | Description |
|------|-------------|
| `h2.xyz` | Hydrogen molecule (nonpolar reference system) |
| `hcl.xyz` | Hydrogen chloride molecule representing a polar covalent bond |
| `h2o.xyz` | Water molecule containing an O–H functional group |
| `nh3.xyz` | Ammonia molecule containing an N–H functional group |
| `methanol.xyz` | Alcohol functional group model (O–H) |
| `formaldehyde.xyz` | Carbonyl functional group model (C=O) |
| `acetaldehyde.xyz` | Aldehyde functional group model |

These molecules were selected to represent chemically and pharmaceutically relevant functional groups used for qualitative comparison of vibrational and polarity trends.

## Final CNDO/2 Results Table
## Final CNDO/2 Results Table

| Molecule | Property | CNDO/2 Result | Experimental / Reference Value | Interpretation |
|---|---|---:|---:|---|
| H₂ | Stretch Vibrational Frequency | 4452 cm⁻¹ | ~4400 cm⁻¹ | Excellent agreement with experimental stretching frequency |
| HCl | Stretch Vibrational Frequency | 2731 cm⁻¹ | ~2991 cm⁻¹ | Reasonable qualitative agreement for polar diatomic system |
| H₂O | Vibrational Mode 1 | 4284 cm⁻¹ | ~1595 cm⁻¹ | Bending mode overestimated by semiempirical approximation |
| H₂O | Vibrational Mode 2 | 4869 cm⁻¹ | ~3657 cm⁻¹ | O–H stretching mode reproduced qualitatively |
| H₂O | Vibrational Mode 3 | 4885 cm⁻¹ | ~3756 cm⁻¹ | High-frequency O–H stretching behavior captured |
| H₂ | Dipole Moment | 0.00 D | 0.00 D | Correctly predicts nonpolar molecular behavior |
| HCl | Dipole Moment | 0.24 D | 1.08 D | Correctly captures molecular polarity trend |
| H₂O | Dipole Moment | 0.45 D | 1.85 D | Polar molecular character reproduced qualitatively |
| NH₃ | Dipole Moment | 0.54 D | 1.47 D | Correctly predicts polarity of amine-containing system |
| Methanol | Dipole Moment | 4.20 D | ~1.70 D | Strong polarity associated with alcohol functional group captured |
| Formaldehyde | Dipole Moment | 2.34 D | ~2.33 D | Excellent agreement for carbonyl polarity |
| Acetaldehyde | Dipole Moment | 8.20 D | ~2.70 D | Polarity trend captured but overestimated for larger organic system |

The CNDO/2 implementation successfully reproduced qualitative molecular trends in both vibrational frequencies and dipole moments. Small diatomic systems such as H₂ showed strong agreement with experimental vibrational frequencies, validating the finite-difference Hessian and mass-weighted vibrational workflow. More complex molecules such as H₂O, NH₃, and methanol exhibited larger deviations from experiment, which is expected for a simplified semiempirical model without geometry optimization or anharmonic corrections. Despite these quantitative limitations, the method correctly captured relative polarity trends and characteristic high-frequency stretching modes, demonstrating that CNDO/2 can provide useful qualitative insight into molecular behavior and pharmaceutically relevant functional-group properties.

The CNDO/2 implementation successfully reproduced qualitative molecular trends in both vibrational frequencies and dipole moments. Small diatomic systems such as H₂ showed strong agreement with experimental vibrational frequencies, validating the finite-difference Hessian and mass-weighted vibrational workflow. More complex molecules such as H₂O, NH₃, and methanol exhibited larger deviations from experiment, which is expected for a simplified semiempirical model without geometry optimization or anharmonic corrections. Despite these quantitative limitations, the method correctly captured relative polarity trends and characteristic high-frequency stretching modes, demonstrating that CNDO/2 can provide useful qualitative insight into molecular behavior and functional-group properties relevant to pharmaceutical chemistry.

## Experimental Validation References (NIST)

Main CCCBDB Database

- NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB)
- https://cccbdb.nist.gov/?utm_source

## Molecule Validation Links

1. H₂ (Hydrogen)
- Experimental Data:
- https://cccbdb.nist.gov/exp2x.asp?casno=1333740&utm

2. HCl (Hydrogen Chloride)
- Experimental Vibrational Data:
- https://cccbdb.nist.gov/exp2x.asp?casno=7647010&charge=0&utm
- Detailed Spectroscopic Constants:
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647010&Mask=20&utm

3. H₂O (Water)
- Experimental Vibrational Frequencies:
- https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0&utm
- Experimental Geometry Data:
- https://cccbdb.nist.gov/expgeom2x.asp?casno=7732185&utm

4. NH₃ (Ammonia)
- Experimental Vibrational Frequencies:
- https://cccbdb.nist.gov/exp2.asp?casno=7664417&utm
- Experimental Geometry Data:
- https://cccbdb.nist.gov/expgeom2x.asp?casno=7664417&charge=0&utm

5. Methanol (CH₃OH)
- Experimental Data:
- https://cccbdb.nist.gov/exp2x.asp?casno=67561&charge=0&utm
- Experimental Geometry:
- https://cccbdb.nist.gov/expgeom2x.asp?casno=67561&charge=0&utm
- NIST WebBook:
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=1000&utm

6. Formaldehyde (CH₂O)
- Experimental Data:
- https://cccbdb.nist.gov/exp2x.asp?casno=50000&charge=0&utm
- Experimental Geometry:
- https://cccbdb.nist.gov/expgeom2x.asp?casno=50000&charge=0&utm
- NIST WebBook:
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C50000&Mask=1000&utm

7. Acetaldehyde (CH₃CHO)
- Experimental Data:
- https://cccbdb.nist.gov/exp2x.asp?casno=75070&charge=0&utm
- Experimental Geometry:
- https://cccbdb.nist.gov/expgeom2x.asp?casno=75070&charge=0&utm
- NIST WebBook:
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070&Mask=1000&utm

## Dipole Moment Validation

Experimental dipole moment reference values used in this project were obtained from the NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB).

NIST Experimental Dipole Moment Database:
https://cccbdb.nist.gov/diplistx.asp

Validated molecules in this project include:
- H2
- HCl
- H2O
- NH3
- Methanol
- Formaldehyde
- Acetaldehyde


Experimental vibrational frequencies, dipole moments, and molecular geometries used for validation in this project were obtained from the National Institute of Standards and Technology (NIST) Computational Chemistry Comparison and Benchmark Database (CCCBDB).
  




## Functional Group Trends

Observed trends include:

- O–H groups produce strong dipole moments and high stretching frequencies
- Carbonyl groups exhibit strong IR-active vibrational modes
- Nonpolar molecules such as H₂ produce negligible dipole moments
- Polar bonds such as H–Cl generate significant molecular polarity

These trends are relevant for understanding molecular behavior in pharmaceutical and biologically relevant systems.

## Limitations

Current limitations of the implementation include:

- Closed-shell systems only
- CNDO/2 semiempirical approximation
- No geometry optimization
- Sensitivity to finite-difference step size
- Reduced quantitative accuracy for larger molecules
- Vibrational analysis most reliable for small molecular systems

This project focuses on reproducing qualitative molecular trends rather than highly accurate ab initio predictions.

## Future Work

Possible future extensions include:

- geometry optimization
- improved SCF convergence schemes
- larger molecular systems
- Density Functional Theory (DFT) comparisons
- solvent effects
- automated visualization pipelines
- expanded pharmaceutical functional-group analysis



University of California, Berkeley

Spring 2026


# References

1. Pople, J. A.; Segal, G. A. Approximate Self‐Consistent Molecular Orbital Theory. II. Calculations with Complete Neglect of Differential Overlap. J. Chem. Phys. 1965.

2. National Institute of Standards and Technology (NIST) Chemistry WebBook:
https://webbook.nist.gov/chemistry/

3. NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB):
https://cccbdb.nist.gov/
