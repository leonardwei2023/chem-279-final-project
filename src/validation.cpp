
#include "validation.h"
#include <iostream>
#include <cmath>

void Validation::compare_with_reference(
    const std::vector<double>& computed,
    const std::vector<double>& reference
) {
    std::cout << "\nValidation vs reference:\n";

    for (size_t i = 0; i < computed.size(); i++) {
        double error = std::abs(computed[i] - reference[i]);

        std::cout << "Mode " << i + 1
                  << " | Computed: " << computed[i]
                  << " | Reference: " << reference[i]
                  << " | Error: " << error << "\n";
    }
}
