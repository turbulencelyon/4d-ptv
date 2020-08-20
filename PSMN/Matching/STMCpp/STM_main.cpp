/*
 * main entry for STM program
 *
 */

#include <ctime>
#include <iostream>
#include <vector>
#include <map>

#include <sys/stat.h>

#include "STM.h"
#include "STM_hdf5.h"
#include "STM_helpers.h"
#include "CLI11.hpp"

void do_STM(const std::string input_dir, const std::string filename, const unsigned int maxframes, const unsigned int mincameras, const double maxdistance, const double multiplematchesperraymindistance, const unsigned int maxmatchesperray, const struct boundingboxspec &bb, bool save_hdf5, bool save_bin, std::string output_dir) {
    std::string inputfile = input_dir + '/' + filename;
    // Print parameters
    std::cout << "Input file: " << inputfile << std::endl;
    std::cout << "maxframes: " << maxframes << std::endl;
    std::cout << "mincameras: " << mincameras << std::endl;
    std::cout << "maxdistance: " << maxdistance << std::endl;
    std::cout << "nx: " << bb.nx << std::endl;
    std::cout << "ny: " << bb.ny << std::endl;
    std::cout << "nz: " << bb.nz << std::endl;
    std::cout << "bounding box: " << bb.xmin << "," << bb.xmax << "," << bb.ymin << "," << bb.ymax << "," << bb.zmin << "," << bb.zmax << std::endl;

    // Go
    std::clock_t tstart = clock();
    
    init();
            
    std::vector<std::vector<double>> bounds;    // Initialize bounds vector
    std::vector<double> tmp;
    tmp.clear();
    for(int i = 0; i<bb.nx;i++)
    {
        tmp.push_back(bb.xmin + i*(bb.xmax-bb.xmin)/bb.nx);
        //std::cout << tmp[tmp.size()-1] << " ";
    }
    tmp.push_back(bb.xmax);
    
    std::cout << "Size:" << tmp.size() << " ";
    bounds.push_back(tmp);
    tmp.clear();
    for(int i = 0; i<bb.ny;i++)
    {
        tmp.push_back(bb.ymin + i*(bb.ymax-bb.ymin)/bb.ny);
        //std::cout << tmp[tmp.size()-1] << " ";
    }
    tmp.push_back(bb.ymax);
    //std::cout << tmp[tmp.size()-1] << " ";
    bounds.push_back(tmp);
    tmp.clear();
    for(int i = 0; i<bb.nz;i++)
    {
        tmp.push_back(bb.zmin + i*(bb.zmax-bb.zmin)/bb.nz);
        //std::cout << tmp[tmp.size()-1] << " ";
    }
    tmp.push_back(bb.zmax);
    //std::cout << tmp[tmp.size()-1] << " ";
    bounds.push_back(tmp);
    
    std::cout << "# of cells: " << bb.nx*bb.ny*bb.nz << "\n";

    std::string filebasename = filename;
    size_t lastindex = filebasename.find_last_of(".");
    filebasename = filebasename.substr(0, lastindex);
    
    // Create output binary file
    std::ofstream streamout;
    if (save_bin) {
        std::string outputfile = output_dir + "/" + filebasename + "_out_cpp.bin";
        streamout.open(outputfile, std::ios::out | std::ios::binary);
        std::cout << "binary output file: " << outputfile << "\n";
    }

    // Create output HDF5 file
    STM_File stm_file;
    if (save_hdf5) {
        std::string outputh5file = output_dir + "/" + filebasename + "_out_cpp.h5";
        stm_file.open(outputh5file, maxframes, mincameras, maxdistance, maxmatchesperray, bb.nx, bb.ny, bb.nz);
        std::cout << "HDF5 outout file: " << outputh5file << std::endl;
    }
    
    // Open input binary file
    std::ifstream rayfile;
    rayfile.open (inputfile, std::ios::in | std::ios::binary);
    
    
    unsigned int currentframe = 0;
    if (rayfile.is_open()) {
        uint32_t numrays = 0;
        rayfile.read((char*)&numrays,sizeof(numrays));
        currentframe++;
        
        while(!rayfile.eof() && currentframe <= maxframes)
        {
            std::cout << "#######\n";
            std::cout << "Frame: " << currentframe << "\nNumber of rays: " << numrays << "\n";
            std::cout << timeasstring() << "\n";
            // Read array of rays
            std::vector<ray> rays;
            for(unsigned int i = 0; i < numrays; i++) {
                ray tmpray;
                uint8_t tmpcamid;
                uint16_t tmprayid;
                float tmpx,tmpy,tmpz,tmpvx,tmpvy,tmpvz;  // Floats (4 byte)
                
                rayfile.read((char*)&tmpcamid,sizeof(tmpcamid));
                rayfile.read((char*)&tmprayid,sizeof(tmprayid));
                rayfile.read((char*)&tmpx,sizeof(tmpx));
                rayfile.read((char*)&tmpy,sizeof(tmpy));
                rayfile.read((char*)&tmpz,sizeof(tmpz));
                rayfile.read((char*)&tmpvx,sizeof(tmpvx));
                rayfile.read((char*)&tmpvy,sizeof(tmpvy));
                rayfile.read((char*)&tmpvz,sizeof(tmpvz));
                tmpray.camid = tmpcamid;
                tmpray.rayid = tmprayid;
                tmpray.x = tmpx;
                tmpray.y = tmpy;
                tmpray.z = tmpz;
                tmpray.vx = tmpvx;
                tmpray.vy = tmpvy;
                tmpray.vz = tmpvz;
                
                rays.push_back(tmpray);
            }
            
            std::vector<candidatematch> results;
    
             // Do the actual stuff, currently no return argument; implement!
            results = SpaceTraversalMatching(rays, bb, bounds, maxmatchesperray, mincameras, maxdistance, multiplematchesperraymindistance);
            
            if (save_bin) {
                // Write results to binary file
                uint32_t numberofmatches = (int)results.size();
                streamout.write((char*)&numberofmatches, sizeof(uint32_t));
                std::cout << "Number of matches: " << numberofmatches << "\n";
                for(auto match: results)
                {
                    uint8_t numberofcams = (int)match.camrayids.size();
                    streamout.write((char*)&numberofcams, sizeof(uint8_t));
                    
                    float val = match.matchx;
                    streamout.write((char*)&val, sizeof(float));
                    val = match.matchy;
                    streamout.write((char*)&val, sizeof(float));
                    val = match.matchz;
                    streamout.write((char*)&val, sizeof(float));
                    val = match.matcherror;
                    streamout.write((char*)&val, sizeof(float));
                    for(auto camrayid:match.camrayids)
                    {
                        uint8_t camid = (int)camrayid.camid;
                        uint16_t rayid =(int)camrayid.rayid;
                        streamout.write((char*)&camid, sizeof(uint8_t));
                        streamout.write((char*)&rayid, sizeof(uint16_t));
                    }
                }
            }

            if (save_hdf5) {
                // Write results to HDF5 file
                stm_file.write_matches(currentframe, results);
            }
            
            rayfile.read((char*)&numrays,sizeof(numrays));  // Read number of rays for next frame
            currentframe++;
        }
        unsigned int lastframe = currentframe - 1;
        if (lastframe != maxframes) {
            std::cout << "Last frame was " << lastframe << " and not " << maxframes << std::endl;
            stm_file.set_last_frame(lastframe);
        }
    }
    else {
        std::cout << "Something went wrong during opening the file. Doez teh file exist?\n";
    }
    rayfile.close();
    streamout.close();
    
    
    std::clock_t tend = clock();
    double timing = double(tend - tstart) / CLOCKS_PER_SEC;
    std::cout << "\nElapsed time = " << timing << " sec\n";
    timing /= currentframe-1;
    std::cout << "\nElapsed time/frame = " << timing << " sec\n";
    
}

int main(int argc, char **argv) {
    CLI::App app{"Stereo-matching program. Created by Sander Huisman."};

    std::string input_file;
    std::string output_dir = ".";
    unsigned int maxframes{0};
    unsigned int mincameras{2};
    unsigned int maxmatchesperray{2};
    double maxdistance{0.2};
    double multiplematchesperraymindistance{0.0};
    unsigned int nx, ny, nz;
    std::vector<int> bounding_box;
    struct boundingboxspec bb;
    bool print_config{false};
    bool hdf5{false};
    bool bin{false};

    app.add_option("-i,--input", input_file, "Input file")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-dir", output_dir, "Output directory");
    app.add_option("-f,--frames", maxframes, "Number of frames")->required(); 
    app.add_option("-c,--mincameras", mincameras, "Minimum number of rays for a match. Default: 2.");
    app.add_option("-d,--maxdistance", maxdistance, "Maximum distance allowed for a match. Default: 0.2.");
    app.add_option("-s,--multiplematchesperraymindistance", multiplematchesperraymindistance, "Minimum allowed distance between matches found for the same ray. Default: 0.0. (allow all).");
    app.add_option("-m,--maxmatchesperray", maxmatchesperray, "Maximum matches per ray. Default: 2.");
    app.add_option("-x,--nx", nx, "nx")->required();
    app.add_option("-y,--ny", ny, "ny")->required();
    app.add_option("-z,--nz", nz, "nz")->required();
    app.add_option("-b,--bb", bounding_box, "Bounding box minx maxx miny maxy minz maxz")->expected(6)->required();
    app.set_config("--config");
    app.add_flag("--print-config", print_config, "Prints an ini config to standard output, and returns.");
    app.add_flag("--hdf5", hdf5, "Save output in HDF5 file");
    app.add_flag("--bin", bin, "Save output in binary file");

    CLI11_PARSE(app, argc, argv);

    if (print_config) {
        std::cout << app.config_to_str(true,true);
        exit(EXIT_SUCCESS);
    }

    if (!hdf5 && !bin) {
        std::cerr << "Select at least one output format, e.g. --hdf5 or --bin." << std::endl;
        exit(EXIT_FAILURE);
    }

    bb.xmin = bounding_box[0];
    bb.xmax = bounding_box[1];
    bb.ymin = bounding_box[2];
    bb.ymax = bounding_box[3];
    bb.zmin = bounding_box[4];
    bb.zmax = bounding_box[5];
    bb.nx = nx;
    bb.ny = ny;
    bb.nz = nz;

    // Check that output_dir does not end with '/'
    if (output_dir.back() == '/') {
        output_dir.pop_back();
    }

    // Check that output_dir exists and create it otherwise
    int status = mkdir(output_dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (status == 0) {
        std::cout << "Directory '" << output_dir << "' created." << std::endl;
    } else if (errno != EEXIST) {
        std::cerr << "Directory '" << output_dir << "' does not exist, and creating it failed with errno=" << errno << std::endl;
    }

    // Split input_file into input_dir and filename
    // 1. check that filename does not end with '/'
    if (input_file.back() == '/') {
        std::cerr << "Input filename should not end with '/'" << std::endl;
        exit(EXIT_FAILURE);
    }
    // 2. look for the last '/' in the string
    std::size_t found = input_file.find_last_of("/");
    std::string input_dir = input_file.substr(0, found);
    std::string filename = input_file.substr(found+1);
    if (found == std::string::npos) {
        input_dir = ".";
        filename = input_file;
    }

    do_STM(input_dir, filename, maxframes, mincameras, maxdistance, multiplematchesperraymindistance, maxmatchesperray, bb, hdf5, bin, output_dir);

    return EXIT_SUCCESS;
}
