/*
 * This file defines the STM_File class which represents a HDF5 file storing the
 * matches found by STM.
 *
 * Header file
 *
 */

#ifndef STM_HDF5_HEADER_H
#define STM_HDF5_HEADER_H

#include <vector>
#include <string>
#include <H5Cpp.h>
#include "STM_types.h"

class STM_File {

    private:

        H5::H5File *f;
        bool readonly_mode;
        bool opened;

    public:

        unsigned int maxframes;
        unsigned int mincameras;
        double maxdistance;
        unsigned int maxmatchesperray;
        unsigned int nx;
        unsigned int ny;
        unsigned int nz;

        // Constructor do nothing
        STM_File();

        // Constructor for STM_File in readonly mode
        STM_File(std::string filename);

        // Constructor for STM_file in write mode
        STM_File(std::string filename, const unsigned int maxframes, const unsigned int mincameras, const double maxdistance, const unsigned int maxmatchesperray, const unsigned int nx, const unsigned int ny, const unsigned int nz);

        // Destructor
        ~STM_File();

        // Instance methods
        void open(std::string filename);
        void open(std::string filename, const unsigned int maxframes, const unsigned int mincameras, const double maxdistance, const unsigned int maxmatchesperray, const unsigned int nx, const unsigned int ny, const unsigned int nz);
        void write_matches(const unsigned int frameno, std::vector<candidatematch> &matches);
        std::vector<candidatematch> read_matches(const unsigned int frameno);
        void set_last_frame(const unsigned int maxframes);
};

#endif
