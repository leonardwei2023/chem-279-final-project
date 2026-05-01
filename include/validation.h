#ifndef VALIDATION_H
#define VALIDATION_H

#include <string>
#include <vector>

class Validation {
public:
    static std::vector<double> read_frequencies(const std::string& filename);

    static void compare(
        const std::vector<double>& computed,
        const std::vector<double>& reference
    );
};

#endif
