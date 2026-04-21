# CNDO/2 Molecular Properties and Vibrational Analysis

CHEM 279 Final project:

Numerical Algorithms Applied to Computational Quantum Chemistry

## Project Overview

In this project, we extend our CNDO/2 implementation from earlier assignments to compute additional molecular properties beyond total energy.

Specifically, we add:

- molecular dipole moment calculation  
- finite-difference Hessian construction  
- vibrational frequency estimation  
- normal mode visualization using `.xyz` animations  

The goal is to show that once we have a converged electronic structure, we can extract meaningful physical properties from it.

## Motivation

Up to now, CNDO/2 has mostly been used to compute energies and densities. However, the density matrix actually contains much more information.

## What this project does

### Dipole moments
We compute the molecular dipole using nuclear + electronic contributions.

### Finite-difference Hessian

H_ij ≈ [E(x+h_i+h_j) - E(x+h_i-h_j) - E(x-h_i+h_j) + E(x-h_i-h_j)] / (4 h^2)

### Vibrational frequencies

We build the Hessian, mass-weight it, and diagonalize.

### Visualization

Normal modes are exported as .xyz animations for Avogadro/VMD.

## Build
```
mkdir build
cd build
cmake ..
cmake --build .
```

## Run

```
./cndo2_project ../examples/hcl.xyz --step 0.01 --animate
```
## Repository Structure
```
chem279_final_project_repo/

├── CMakeLists.txt        
├── README.md
├── .gitignore             

├── src/                  
│   ├── main.cpp          
│   ├── molecule.cpp      
│   ├── dipole.cpp        
│   ├── finite_difference.cpp  
│   ├── modes.cpp         
│   └── xyz_io.cpp        

├── include/              
│   ├── molecule.hpp
│   ├── dipole.hpp
│   ├── finite_difference.hpp
│   ├── modes.hpp
│   ├── xyz_io.hpp
│   └── cndo_engine.hpp   

├── examples/            
│   ├── h2.xyz
│   ├── hcl.xyz
│   ├── h2o.xyz
│   └── nh3.xyz


├── results/              
│   ├── dipoles.txt
│   ├── frequencies.txt
│   └── modes.xyz

└── build/                
```

## Notes

## Authors

David Houshangi

Leonard Ming Wei
