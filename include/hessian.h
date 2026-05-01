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

    void resize(int new_size);

    void read_from_file(const std::string& filename, int expected_size);
    void write_to_file(const std::string& filename) const;

    int get_size() const;

    double get_value(int i, int j) const;
    void set_value(int i, int j, double value);

    std::vector<std::vector<double>> get_matrix() const;
};

#endif
