/*
 * Common STM type definitions
 *
 */

#ifndef STM_HEADER_H
#define STM_HEADER_H

#include <vector>

struct cellvisit{
    double time;
    int dimension = -1;
    bool edgeQ = false;
};

struct cellid{
    long xi;
    long yi;
    long zi;
};

struct traversedcell{
    int camid = -1;
    int rayid = -1;
    struct cellid cellid;
};

struct camrayid{
    long camid;
    long rayid;
};

struct groupedcell{
    struct cellid cellid;
    std::vector<camrayid> camrayids;
};

struct candidatematch{
    double matchx = 0;
    double matchy = 0;
    double matchz = 0;
    double matcherror = -1;
    std::vector<camrayid> camrayids;
};

struct ray{
    int camid = -1;
    int rayid = -1;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

struct transformedray{
    int camid = -1;
    int rayid = -1;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    bool hit = false;
    bool inside = false;
};

struct boundingboxspec{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    int nx;
    int ny;
    int nz;
};

struct hitpoint{
    double t;
    double posx;
    double posy;
    double posz;
};

#endif
