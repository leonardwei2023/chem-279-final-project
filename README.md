# Calculation of Molecular Dipole Moments and Basic Vibrational Frequencies using CNDO/2. 

David Houshangi - davidhoushangi@berkeley.edu 

Leonard Ming Wei - dranoelmi@berkeley.edu  

# Project Overview 

Our goal is to extend our current CNDO/2 density calculations to include dipole moments and to 
visualize the vibrational frequencies of a molecule. 

# Objectives 

1. Our first objective is to compute the molecular dipole. 
2. Our second objective is to implement a basic implementation of vibrational analysis.

# Methodology 

## a. Dipole Calculations 

i. Compute x, y, z dipole components from density matrix and atomic partial 
charges 

ii. Compute final dipole vector from components and convert to Debye

iii. (Optional) Visualize dipole vector and charge distribution using VMD or 
Avogadro 

## b. Vibrational Analysis (Finite Differences) 

i. Randomly displace atoms 

ii. Recompute CNDO/2 and Hessian 

iii. Convert Hessian to mass-weighted coordinates  𝐻𝑖𝑗
' = 𝐻𝑖𝑗
𝑚
𝑖
𝑚𝑗

iv. Diagonalize mass-weighted Hessian to extract frequencies  𝑣 𝑖 
= 1

v. Plot frequencies using matplotlib 


# Planned Experiments 
2π𝑐
λ
𝑖
These implementations will be tested on small molecules with simpler geometric structures to 
ensure numerical accuracy. 

1. Dipole moments will be calculated for molecules such as Hydrogen Chloride (HCl), 
water (H2O), ammonium (NH3).

2. Vibrational frequency analysis will be performed on small molecules such as HCl and H2 
to maintain the computational efficiency.
  
Upon analysis of the results, trends in molecular polarity and vibrational behavior will be 
identified, with attempts for them to be qualitatively juxtaposed to expected physical behavior. 
