#include "hessian.h"

#include <fstream>
#include <stdexcept>

Hessian::Hessian() {
    size = 0;
}

void Hessian::resize(int new_size) {
    size = new_size;
    matrix.clear();
    matrix.resize(size, std::vector<double>(size, 0.0));
}

void Hessian::read_from_file(const std::string& filename, int expected_size) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open Hessian file: " + filename);
    }

    resize(expected_size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!(file >> matrix[i][j])) {
                throw std::runtime_error("Hessian file does not contain enough values.");
            }
        }
    }
}

void Hessian::write_to_file(const std::string& filename) const {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not write Hessian file: " + filename);
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            file << matrix[i][j];

            if (j < size - 1) {
                file << " ";
            }
        }

        file << "\n";
    }
}

int Hessian::get_size() const {
    return size;
}

double Hessian::get_value(int i, int j) const {
    return matrix[i][j];
}

void Hessian::set_value(int i, int j, double value) {
    matrix[i][j] = value;
}

std::vector<std::vector<double>> Hessian::get_matrix() const {
    return matrix;
}
