
#include "cndo_engine.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>

CNDOEngine::CNDOEngine(const std::string& command_in) {
    command = command_in;
    file_counter = 0;
}

std::string CNDOEngine::make_temp_filename() {
    file_counter++;
    return "tmp_cndo_displaced_" + std::to_string(file_counter) + ".xyz";
}

double CNDOEngine::extract_energy_from_output(const std::string& output) const {
    std::stringstream ss(output);
    std::string token;

    bool found_number = false;
    double last_number = 0.0;

    while (ss >> token) {
        char* end_pointer = nullptr;
        double value = std::strtod(token.c_str(), &end_pointer);

        if (end_pointer != token.c_str()) {
            found_number = true;
            last_number = value;
        }
    }

    if (!found_number) {
        throw std::runtime_error("Could not find a numeric energy in CNDO/2 output.");
    }

    return last_number;
}

double CNDOEngine::compute_energy(const Molecule& molecule) {
    std::string temp_file = make_temp_filename();

    molecule.write_xyz(temp_file);

    std::string full_command = command + " " + temp_file;

    FILE* pipe = popen(full_command.c_str(), "r");

    if (pipe == nullptr) {
        throw std::runtime_error("Could not run CNDO/2 energy command.");
    }

    char buffer[256];
    std::string output;

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        output += buffer;
    }

    int status = pclose(pipe);

    std::remove(temp_file.c_str());

    if (status != 0) {
        throw std::runtime_error("CNDO/2 energy command failed.");
    }

    return extract_energy_from_output(output);
}
