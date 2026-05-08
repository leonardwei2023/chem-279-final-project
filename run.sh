
#!/bin/bash

set -e

echo "CHEM 279 Final Project Run Script"
echo "CNDO/2 Vibrational Frequencies and Dipoles"

echo ""
echo "Creating build directory..."
mkdir -p build
cd build

echo ""
echo "Configuring project with CMake..."
cmake ..

echo ""
echo "Building project..."
cmake --build .

echo ""
echo "Setting CNDO/2 energy executable..."
export CNDO_ENERGY_CMD="./cndo_energy"

echo ""
echo "Running vibrational frequency calculations..."

echo "H2 vibration..."
./vibrational_frequency finite-diff-vib ../input/h2.xyz h2_fd.dat 0.005

echo "HCl vibration..."
./vibrational_frequency finite-diff-vib ../input/hcl.xyz hcl_fd.dat 0.005

echo "H2O vibration..."
./vibrational_frequency finite-diff-vib ../input/h2o.xyz h2o_fd.dat 0.005

echo ""
echo "Running dipole moment calculations..."

echo "H2 dipole..."
./vibrational_frequency dipole ../input/h2.xyz

echo "HCl dipole..."
./vibrational_frequency dipole ../input/hcl.xyz

echo "H2O dipole..."
./vibrational_frequency dipole ../input/h2o.xyz

echo "NH3 dipole..."
./vibrational_frequency dipole ../input/nh3.xyz

echo "Methanol dipole..."
./vibrational_frequency dipole ../input/methanol.xyz

echo "Formaldehyde dipole..."
./vibrational_frequency dipole ../input/formaldehyde.xyz

echo "Acetaldehyde dipole..."
./vibrational_frequency dipole ../input/acetaldehyde.xyz

echo ""
echo "Generating optional normal mode animation for H2..."
./vibrational_frequency finite-diff-vib ../input/h2.xyz h2_fd.dat 0.005 --animate

echo ""
echo " Chem 279 Final Project Report - Done."
echo "All calculations completed successfully."
echo "Generated files are located in the build directory."
echo "Examples:"
echo "  h2_fd.dat"
echo "  hcl_fd.dat"
echo "  h2o_fd.dat"
echo "  computed_frequencies.dat"
echo "  dipole_moment.dat"
echo "  normal_modes.xyz"
echo "  Thank you Dr. Agrawal and Lizzie."
echo "  End."
