/*
 * This file defines the STM_File class which represents a HDF5 file storing the
 * matches found by STM.
 *
 * Implementation file
 *
 */

#include <iostream>
#include "STM_hdf5.h"

STM_File::STM_File() {
    opened = false;
}

STM_File::STM_File(std::string filename) {
    opened = false;
    open(filename);
}

STM_File::STM_File(std::string filename, const unsigned int maxframes, const unsigned int mincameras, const double maxdistance, const unsigned int maxmatchesperray, const unsigned int nx, const unsigned int ny, const unsigned int nz) {
    opened = false;
    open(filename, maxframes, mincameras, maxdistance, maxmatchesperray, nx, ny, nz);
}

void STM_File::open(std::string filename) {
    /*
     * Constructor method
     * Open STM_file in readonly mode
     *
     */

    if (opened) {
        std::cerr << "Warning: STM HDF5 file was already open!" << std::endl;
    }

    f = new H5::H5File(H5std_string(filename), H5F_ACC_RDONLY);
    readonly_mode = true;
    std::cout << "File '" << filename << "' opened in readonly mode." << std::endl;

    // Read the attributes from the HDF5 file
    H5::Attribute attr_maxframes = f->openAttribute(H5std_string("maxframes"));
    attr_maxframes.read(H5::PredType::NATIVE_INT, &maxframes);
    std::cout << "maxframes = " << maxframes << std::endl;

    H5::Attribute attr_mincameras = f->openAttribute(H5std_string("mincameras"));
    attr_mincameras.read(H5::PredType::NATIVE_INT, &mincameras);
    std::cout << "mincameras = " << mincameras << std::endl;

    H5::Attribute attr_maxdistance = f->openAttribute(H5std_string("maxdistance"));
    attr_maxdistance.read(H5::PredType::NATIVE_DOUBLE, &maxdistance);
    std::cout << "maxdistance = " << maxdistance << std::endl;

    H5::Attribute attr_maxmatchesperray = f->openAttribute(H5std_string("maxmatchesperray"));
    attr_maxmatchesperray.read(H5::PredType::NATIVE_INT, &maxmatchesperray);
    std::cout << "maxmatchesperray = " << maxmatchesperray << std::endl;

    H5::Attribute attr_nx = f->openAttribute(H5std_string("nx"));
    attr_nx.read(H5::PredType::NATIVE_INT, &nx);
    std::cout << "nx = " << nx << std::endl;

    H5::Attribute attr_ny = f->openAttribute(H5std_string("ny"));
    attr_ny.read(H5::PredType::NATIVE_INT, &ny);
    std::cout << "ny = " << ny << std::endl;

    H5::Attribute attr_nz = f->openAttribute(H5std_string("nz"));
    attr_nz.read(H5::PredType::NATIVE_INT, &nz);
    std::cout << "nz = " << nz << std::endl;

    opened = true;
}

void STM_File::open(std::string filename, const unsigned int maxframes, const unsigned int mincameras, const double maxdistance, const unsigned int maxmatchesperray, const unsigned int nx, const unsigned int ny, const unsigned int nz) {
    /*
     * Constructor method
     * Open STM_file in write mode
     *
     */

    if (opened) {
        std::cerr << "Warning: STM HDF5 file was already open!" << std::endl;
    }

    f = new H5::H5File(H5std_string(filename), H5F_ACC_TRUNC);
    readonly_mode = false;
    std::cout << "File '" << filename << "' opened in write mode." << std::endl;

    this->maxframes = maxframes;
    this->mincameras = mincameras;
    this->maxdistance = maxdistance;
    this->maxmatchesperray = maxmatchesperray;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    // Write to HDF5 file
    const unsigned int RANK = 1;
    hsize_t dim_vector[RANK] = {1};
    H5::DataSpace attr_dataspace = H5::DataSpace(RANK, dim_vector);

    H5::Attribute attr_maxframes = f->createAttribute(H5std_string("maxframes"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_maxframes.write(H5::PredType::NATIVE_INT, &maxframes);

    H5::Attribute attr_mincameras = f->createAttribute(H5std_string("mincameras"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_mincameras.write(H5::PredType::NATIVE_INT,&mincameras);

    H5::Attribute attr_maxdistance = f->createAttribute(H5std_string("maxdistance"),
                                                      H5::PredType::NATIVE_DOUBLE,
                                                      attr_dataspace);
    attr_maxdistance.write(H5::PredType::NATIVE_DOUBLE, &maxdistance);

    H5::Attribute attr_maxmatchesperray = f->createAttribute(H5std_string("maxmatchesperray"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_maxmatchesperray.write(H5::PredType::NATIVE_INT, &maxmatchesperray);

    H5::Attribute attr_nx = f->createAttribute(H5std_string("nx"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_nx.write(H5::PredType::NATIVE_INT, &nx);

    H5::Attribute attr_ny = f->createAttribute(H5std_string("ny"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_ny.write(H5::PredType::NATIVE_INT, &ny);

    H5::Attribute attr_nz = f->createAttribute(H5std_string("nz"),
                                                      H5::PredType::NATIVE_INT,
                                                      attr_dataspace);
    attr_nz.write(H5::PredType::NATIVE_INT, &nz);

    opened = true;
}

// Destructor method
STM_File::~STM_File() {
    delete f;
}

// Methods

void STM_File::write_matches(const unsigned int frameno, std::vector<candidatematch> &matches) {
    /*
     * Write matches to HDF5 file.
     *
     */

    // Find the largest number of camrayids for this frame
    unsigned int maxcams = 0;
    for(auto match: matches) {
        if (match.camrayids.size() > maxcams) {
            maxcams = match.camrayids.size();
        }
    }

    // Dataset sizes
    const unsigned int RANK = 2;
    hsize_t xyze_dim_vector[RANK] = {4, matches.size()};
    H5::DataSpace xyze_space(RANK, xyze_dim_vector);
    hsize_t camrayids_dim_vector[RANK] = {2*maxcams, matches.size()};
    H5::DataSpace camrayids_space(RANK, camrayids_dim_vector);

    // Dataset options
    unsigned int chunksize = 1024;
    if (chunksize >= matches.size()) {
        chunksize = matches.size()/4;
    }
    H5::DSetCreatPropList xyze_plist;
    hsize_t xyze_chunksize[RANK] = {4, chunksize};
    if (chunksize > 1) {
        xyze_plist.setChunk(RANK, xyze_chunksize);
        xyze_plist.setDeflate(9);
    }

    H5::DSetCreatPropList camrayids_plist;
    hsize_t camrayids_chunksize[RANK] = {2*maxcams, chunksize};
    if (chunksize > 1) {
        camrayids_plist.setChunk(RANK, camrayids_chunksize);
        camrayids_plist.setDeflate(9);
    }

    // Create the datasets
    std::string xyze_name = "/frame";
    xyze_name.append(std::to_string(frameno));
    xyze_name.append("_xyze");
    H5::DataSet xyze_dset = f->createDataSet(H5std_string(xyze_name),
                                            H5::PredType::NATIVE_DOUBLE,
                                            xyze_space, xyze_plist);
    std::string camrayids_name = "/frame";
    camrayids_name.append(std::to_string(frameno));
    camrayids_name.append("_camrayids");
    H5::DataSet camrayids_dset = f->createDataSet(H5std_string(camrayids_name),
                                                  H5::PredType::NATIVE_LONG,
                                                  camrayids_space, camrayids_plist);

    // Write the data
    hsize_t row_xyze_dim_vector[RANK] = {4, 1};
    H5::DataSpace row_xyze_space(RANK, row_xyze_dim_vector);
    hsize_t row_camrayids_dim_vector[RANK] = {2*maxcams, 1};
    H5::DataSpace row_camrayids_space(RANK, row_camrayids_dim_vector);

    for(unsigned int row = 0; row < matches.size(); row++) {
        // Prepare row data
        double xyze_data[4] = {matches[row].matchx, 
                               matches[row].matchy, 
                               matches[row].matchz, 
                               matches[row].matcherror};
        long camrayids_data[2*maxcams];
        for(unsigned int idno = 0; idno < maxcams; idno++) {
            if (matches[row].camrayids.size() > idno) {
                camrayids_data[2*idno] = matches[row].camrayids[idno].camid;
                camrayids_data[2*idno+1] = matches[row].camrayids[idno].rayid;
            } else {
                camrayids_data[2*idno] = -1;
                camrayids_data[2*idno+1] = -1;
            }
        }

        // Write row
        hsize_t offset[2] = {0, row};

        H5::DataSpace xyze_subspace = xyze_dset.getSpace();
        xyze_subspace.selectHyperslab(H5S_SELECT_SET, row_xyze_dim_vector, offset);
        xyze_dset.write(xyze_data, H5::PredType::NATIVE_DOUBLE, row_xyze_space, xyze_subspace);

        H5::DataSpace camrayids_subspace = camrayids_dset.getSpace();
        camrayids_subspace.selectHyperslab(H5S_SELECT_SET, row_camrayids_dim_vector, offset);
        camrayids_dset.write(camrayids_data, H5::PredType::NATIVE_LONG, row_camrayids_space, camrayids_subspace);

    }
}

std::vector<candidatematch> STM_File::read_matches(const unsigned int frameno) {
    // Open datasets
    std::string xyze_name = "/frame";
    xyze_name.append(std::to_string(frameno));
    xyze_name.append("_xyze");
    H5::DataSet xyze_dset = f->openDataSet(H5std_string(xyze_name));

    std::string camrayids_name = "/frame";
    camrayids_name.append(std::to_string(frameno));
    camrayids_name.append("_camrayids");
    H5::DataSet camrayids_dset = f->openDataSet(H5std_string(camrayids_name));

    // Dataset sizes
    const unsigned int RANK = 2;
    H5::DataSpace xyze_space = xyze_dset.getSpace();
    if (xyze_space.getSimpleExtentNdims() != RANK) {
        std::cerr << "Wrong RANK in dataset!" << std::endl;
        exit(EXIT_FAILURE);
    }
    hsize_t xyze_dim_vector[RANK];
    xyze_space.getSimpleExtentDims(xyze_dim_vector, NULL);
    if (xyze_dim_vector[0] != 4) {
        std::cerr << "Wrong xyze dataset dimension" << std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout << "number of rows in xyze = " << xyze_dim_vector[1] << std::endl;

    H5::DataSpace camrayids_space = camrayids_dset.getSpace();
    if (camrayids_space.getSimpleExtentNdims() != RANK) {
        std::cerr << "Wrong RANK in dataset!" << std::endl;
        exit(EXIT_FAILURE);
    }
    hsize_t camrayids_dim_vector[RANK];
    camrayids_space.getSimpleExtentDims(camrayids_dim_vector, NULL);
    // std::cout << "number of rows in camrayids = " << camrayids_dim_vector[1] << std::endl;
    unsigned int maxcams = camrayids_dim_vector[0]/2;
    if (xyze_dim_vector[1] != camrayids_dim_vector[1]) {
        std::cerr << "xyze camrayids dimension mismatch" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read data
    const unsigned int nrows = xyze_dim_vector[1];
    std::vector<candidatematch> matches(nrows);
    hsize_t row_xyze_dim_vector[RANK] = {4, 1};
    H5::DataSpace row_xyze_space(RANK, row_xyze_dim_vector);
    hsize_t row_camrayids_dim_vector[RANK] = {2*maxcams, 1};
    H5::DataSpace row_camrayids_space(RANK, row_camrayids_dim_vector);

    for(unsigned int row = 0; row < nrows; row++) {
        double xyze_data[4];
        long camrayids_data[2*maxcams];

        hsize_t offset[2] = {0, row};
        H5::DataSpace xyze_subspace = xyze_dset.getSpace();
        xyze_subspace.selectHyperslab(H5S_SELECT_SET, row_xyze_dim_vector, offset);
        xyze_dset.read(xyze_data, H5::PredType::NATIVE_DOUBLE, row_xyze_space, xyze_subspace);

        H5::DataSpace camrayids_subspace = camrayids_dset.getSpace();
        camrayids_subspace.selectHyperslab(H5S_SELECT_SET, row_camrayids_dim_vector, offset);
        camrayids_dset.read(camrayids_data, H5::PredType::NATIVE_LONG, row_camrayids_space, camrayids_subspace);
        matches[row].matchx = xyze_data[0];
        matches[row].matchy = xyze_data[1];
        matches[row].matchz = xyze_data[2];
        matches[row].matcherror = xyze_data[3];

        for(unsigned int idno = 0; idno < maxcams; idno++) {
            struct camrayid camrayid;
            camrayid.camid = camrayids_data[2*idno];
            camrayid.rayid = camrayids_data[2*idno+1];
            if (camrayid.camid != -1 && camrayid.rayid != -1) {
                matches[row].camrayids.push_back(camrayid);
            }
        }
    }

    return matches;
}

void STM_File::set_last_frame(const unsigned int maxframes) {
    /*
     * Override maxframes attribute
     *
     */

    this->maxframes = maxframes;
    if (opened && !readonly_mode) {
        H5::Attribute attr_maxframes = f->openAttribute(H5std_string("maxframes"));
        attr_maxframes.write(H5::PredType::NATIVE_INT, &maxframes);
    }
}
