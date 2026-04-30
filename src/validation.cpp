#include "validation.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

std::vector<double> Validation::read_frequencies(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open frequency file: " + filename);
    }

    std::vector<double> frequencies;
    double value;

    while (file >> value) {
        frequencies.push_back(value);
    }

    return frequencies;
}

void Validation::compare(
    const std::vector<double>& computed,
    const std::vector<double>& reference
) {
    std::cout << "\nValidation against Psi4/reference:\n";

    size_t n = computed.size();

    if (reference.size() < n) {
        n = reference.size();
    }

    for (size_t i = 0; i < n; i++) {
        double error = computed[i] - reference[i];
        double abs_error = std::abs(error);

        std::cout << "Mode " << i + 1
                  << " | Computed: " << computed[i]
                  << " cm^-1 | Reference: " << reference[i]
                  << " cm^-1 | Error: " << abs_error
                  << " cm^-1\n";
    }
}
