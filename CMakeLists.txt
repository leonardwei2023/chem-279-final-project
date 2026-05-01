cmake_minimum_required(VERSION 3.10)
project(cndo2_vibrational)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 (header-only linear algebra library)
find_package(Eigen3 REQUIRED)

# Include project headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# Main executable
# vibrational_frequency: vibration modes, dipole moment, validation
add_executable(vibrational_frequency
    src/main.cpp
    src/molecule.cpp
    src/hessian.cpp
    src/cndo_engine.cpp
    src/finite_difference.cpp
    src/vibration.cpp
    src/validation.cpp
    src/dipole.cpp
)
target_link_libraries(vibrational_frequency Eigen3::Eigen)

# CNDO/2 energy oracle 
# cndo_energy: standalone binary called by CNDOEngine for finite-difference Hessian
add_executable(cndo_energy
    src/cndo_energy.cpp
    src/molecule.cpp
)
target_link_libraries(cndo_energy Eigen3::Eigen)
