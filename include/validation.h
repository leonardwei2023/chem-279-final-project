
#ifndef VALIDATION_H
#define VALIDATION_H

#include <vector>

class Validation {
public:
    static void compare_with_reference(
        const std::vector<double>& computed,
        const std::vector<double>& reference
    );
};

#endif
