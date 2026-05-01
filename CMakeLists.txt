cmake_minimum_required(VERSION 3.10)

project(vibrational_frequency)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Include headers
include_directories(include)

# Main executable (ONLY ONE — no cndo_energy anymore)
add_executable(vibrational_frequency
    src/main.cpp
    src/molecule.cpp
    src/hessian.cpp
    src/cndo_engine.cpp
    src/finite_difference.cpp
    src/vibration.cpp
    src/validation.cpp
)

# Link libraries
target_link_libraries(vibrational_frequency Eigen3::Eigen)
