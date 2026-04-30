cmake_minimum_required(VERSION 3.10)

project(vibrational_frequency)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

include_directories(include)

add_executable(vibrational_frequency
    src/main.cpp
    src/molecule.cpp
    src/hessian.cpp
    src/vibration.cpp
)

target_link_libraries(vibrational_frequency Eigen3::Eigen)
