# Vibrational Frequencies and Dipole Moment using CNDO/2

## Overview

This project implements a computational chemistry workflow in C++ to calculate:

- **Vibrational frequencies** using the **finite difference method**
- **Dipole moments** using a **CNDO/2 self-consistent field (SCF) method**

The purpose of this project is to demonstrate how molecular properties can be computed using numerical techniques and semi-empirical quantum chemistry models.

---

## Features

- CNDO/2 SCF energy calculation
- Finite-difference Hessian construction
- Vibrational frequency computation via eigenvalue analysis
- Dipole moment calculation from electronic structure
- Modular C++ design
- Linear algebra using Eigen

---

## Project Structure
```
vibrational_frequency/
│
├── include/
│ ├── molecule.h
│ ├── hessian.h
│ ├── cndo_engine.h
│ ├── finite_difference.h
│ ├── vibration.h
│ └── validation.h
│
├── src/
│ ├── main.cpp
│ ├── molecule.cpp
│ ├── hessian.cpp
│ ├── cndo_engine.cpp
│ ├── finite_difference.cpp
│ ├── vibration.cpp
│ └── validation.cpp
│
├── CMakeLists.txt
└── README.md

```


---

## Methods

### 1. CNDO/2 SCF Energy

The electronic energy is computed using the CNDO/2 approximation:

- Construct Fock matrix
- Solve for molecular orbitals
- Update density matrix
- Iterate until convergence (SCF cycle)

---

### 2. Finite Difference Hessian

The Hessian matrix is computed numerically using central finite differences:
H_ij ≈ [E(x_i + h, x_j + h) - E(x_i + h, x_j - h)
- E(x_i - h, x_j + h) + E(x_i - h, x_j - h)] / (4h^2)


---

### 3. Vibrational Frequencies

Steps:

1. Construct mass-weighted Hessian
2. Compute eigenvalues
3. Convert to frequencies:

ω = sqrt(λ)


where λ are eigenvalues of the mass-weighted Hessian.

---

### 4. Dipole Moment

The dipole moment is computed as:


μ = Σ (q_i * r_i)


where:
- q_i = atomic charge
- r_i = atomic position

Charges are obtained from the CNDO/2 density matrix.

---

## Build Instructions

### Requirements

- C++17 compiler
- CMake ≥ 3.10
- Eigen library

### Build

```bash
mkdir build
cd build
cmake ..
make
```

# Run
```
./vibrational_frequency
```

# Output

## The program outputs:

- Total SCF energy
- Dipole moment (x, y, z components)
- Vibrational frequencies

# Validation

## Results can be validated using external quantum chemistry packages such as:

- Psi4
- Gaussian

# Notes
- This is a simplified educational implementation
- CNDO/2 is an approximate method
- Accuracy depends on finite difference step size

  
# Future Work
- Add analytical gradients
- Improve SCF convergence
- Add molecular input file support
- Visualize vibrational modes
- Parallelize computations

# Author
