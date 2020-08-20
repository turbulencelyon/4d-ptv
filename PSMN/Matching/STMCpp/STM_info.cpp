/*
 * This program opens an existing STM HDF5 file and prints information.
 *
 */

#include <vector>
#include <iostream>

#include "STM_types.h"
#include "STM_hdf5.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "Usage: STM_info filename" << std::endl;
        exit(EXIT_FAILURE);
    }
    STM_File stm(argv[1]);
    for (unsigned int frameno = 1; frameno <= stm.maxframes; frameno++) {
        std::vector<candidatematch> matches = stm.read_matches(frameno);
        std::cout << "Frame " << frameno << std::endl;
        std::cout << "    No matches: " << matches.size() << std::endl;
    }
}

