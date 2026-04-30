
#include "hessian.h"

#include <fstream>
#include <stdexcept>

Hessian::Hessian() {
    size = 0;
}

void Hessian::read_from_file(const std::string& filename, int expected_size) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open Hessian file: " + filename);
    }

    size = expected_size;
    matrix.clear();
    matrix.resize(size, std::vector<double>(size, 0.0));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!(file >> matrix[i][j])) {
                throw std::runtime_error("Hessian file does not have enough values.");
            }
        }
    }

    file.close();
}

int Hessian::get_size() const {
    return size;
}

std::vector<std::vector<double>> Hessian::get_matrix() const {
    return matrix;
}
