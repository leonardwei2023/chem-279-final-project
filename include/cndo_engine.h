#ifndef CNDO_ENGINE_H
#define CNDO_ENGINE_H

#include "molecule.h"

#include <string>

class CNDOEngine {
private:
    std::string command;
    int file_counter;

    std::string make_temp_filename();
    double extract_energy_from_output(const std::string& output) const;

public:
    CNDOEngine(const std::string& command_in);

    double compute_energy(const Molecule& molecule);
};

#endif
