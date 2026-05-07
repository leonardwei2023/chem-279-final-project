// validation.h
#ifndef VALIDATION_H
#define VALIDATION_H

#include <string>
#include <vector>

class Validation {
public:
    // Reads a plain text file containing one number per line.
    // This is used for both frequency files and p_diagonal.dat.
    static std::vector<double> read_values(const std::string& filename);

    // Compare computed vibrational frequencies with reference values.
    static void compare(
        const std::vector<double>& computed,
        const std::vector<double>& reference
    );
};

#endif
