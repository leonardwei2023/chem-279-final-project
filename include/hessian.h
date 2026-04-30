
#ifndef HESSIAN_H
#define HESSIAN_H

#include <string>
#include <vector>

class Hessian {
private:
    int size;
    std::vector<std::vector<double>> matrix;

public:
    Hessian();

    void read_from_file(const std::string& filename, int expected_size);

    int get_size() const;
    std::vector<std::vector<double>> get_matrix() const;
};

#endif
