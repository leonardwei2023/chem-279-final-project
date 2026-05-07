#include "validation.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> Validation::read_values(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<double> values;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        try {
            values.push_back(std::stod(line));
        } catch (...) {
            // Ignore non-numeric lines.
        }
    }

    if (values.empty()) {
        throw std::runtime_error("No numeric values found in file: " + filename);
    }

    return values;
}

void Validation::compare(
    const std::vector<double>& computed,
    const std::vector<double>& reference
) {
    if (computed.empty()) {
        throw std::runtime_error("Computed frequency list is empty.");
    }

    if (reference.empty()) {
        throw std::runtime_error("Reference frequency list is empty.");
    }

    size_t n = std::min(computed.size(), reference.size());

    double total_abs_error = 0.0;
    double total_percent_error = 0.0;

    std::cout << "\n--- Validation against reference frequencies ---\n";

    for (size_t i = 0; i < n; i++) {
        double error = computed[i] - reference[i];
        double abs_error = std::abs(error);

        double percent_error = 0.0;
        if (std::abs(reference[i]) > 1.0e-12) {
            percent_error = 100.0 * abs_error / std::abs(reference[i]);
        }

        total_abs_error += abs_error;
        total_percent_error += percent_error;

        std::cout << "Mode " << i + 1
                  << " | Computed: " << computed[i]
                  << " cm^-1 | Reference: " << reference[i]
                  << " cm^-1 | Abs Error: " << abs_error
                  << " cm^-1 | Percent Error: " << percent_error
                  << "%\n";
    }

    std::cout << "-----------------------------------------------\n";
    std::cout << "Compared modes: " << n << "\n";
    std::cout << "Mean absolute error: "
              << total_abs_error / static_cast<double>(n)
              << " cm^-1\n";
    std::cout << "Mean percent error: "
              << total_percent_error / static_cast<double>(n)
              << "%\n";

    if (computed.size() != reference.size()) {
        std::cout << "Note: computed and reference files contain different "
                  << "numbers of modes.\n";
    }
}
