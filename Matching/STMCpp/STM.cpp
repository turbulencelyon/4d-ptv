//  main.cpp
//  STM
//
//  Created by Sander Huisman on 10/05/2017.

#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <tuple>

#include "STM.h"
#include "STM_CMatrix.h"
#include "STM_hdf5.h"
#include "STM_helpers.h"

// Global variables

int firstround = 0;
long calls = 0;
std::vector<std::vector<int>> neighbours;   // Stores the relative indices of neighbours
std::vector<std::vector<double>> bounds;

// Implementation

bool comparecellvisitbytime(const cellvisit &a, const cellvisit &b)
{
    return a.time < b.time;
}

bool comparecandidatematches(const candidatematch &a, const candidatematch &b) // First sort by number of cams then by match error
{
    return(a.camrayids.size() > b.camrayids.size() || (a.camrayids.size() == b.camrayids.size() && a.matcherror < b.matcherror));
}

bool comparecamrayidscam(const camrayid &a, const camrayid &b)
{
    return(a.camid < b.camid);
}

bool comparecamrayididentical(const camrayid &a, const camrayid &b)
{
    return(a.camid == b.camid && a.rayid == b.rayid);
}

bool comparecamrayidscamray(const camrayid &a, const camrayid &b)
{
    return(a.camid < b.camid || (a.camid == b.camid && a.rayid < b.rayid));
}

bool comparecamrayidsamecam(const camrayid &a, const camrayid &b)
{
    return(a.camid == b.camid);
}

bool comparecamrayidsidentical(const std::vector<camrayid> &a, const std::vector<camrayid> &b)
{
    if(a.size() != b.size())
    {
        return(false);
    }
    else
    {
        for(unsigned int i = 0; i < a.size();i++)
        {
            if(a[i].camid != b[i].camid || a[i].rayid != b[i].rayid)
            {
                return(false);
            }
        }
        return(true);
    }
}

bool comparecamrayidsordered(const std::vector<camrayid> &a, const std::vector<camrayid> &b)
{
    size_t as = a.size();
    size_t bs = b.size();
    size_t minlen = std::min(as,bs);
    for(size_t i = 0; i < minlen; i++)
    {
        if(a[i].camid < b[i].camid || (a[i].camid == b[i].camid && a[i].rayid < b[i].rayid))
        {
            return(true);
        }
        else
        {
            if(a[i].camid > b[i].camid || (a[i].camid == b[i].camid && a[i].rayid > b[i].rayid))
            {
                return(false);
            }
        }
    }
    if(as != bs)
    {
        return(as < bs);    // Short goes before long
    }
    else                    // This basically means that they have same length and content
    {
        return(false);
    }
}

bool comparetraversedcell(const traversedcell &a, const traversedcell &b)
{
    return(a.cellid.xi < b.cellid.xi || (a.cellid.xi == b.cellid.xi && a.cellid.yi < b.cellid.yi) || (a.cellid.xi == b.cellid.xi && a.cellid.yi == b.cellid.yi && a.cellid.zi < b.cellid.zi));
}

bool comparetraversedcellsamecell(const traversedcell &a, const traversedcell &b)
{
    return(a.cellid.xi == b.cellid.xi && a.cellid.yi == b.cellid.yi && a.cellid.zi == b.cellid.zi);
}

bool comparehitpointbytime(const hitpoint &a, const hitpoint &b)
{
    return a.t < b.t;
}

long positionsorted(std::vector<double> boundaries,double n)
{
    long mn = 0;
    long mx = boundaries.size() - 1;
    long trial;
    boundaries[mn] -= 1.0e-8;
    boundaries[mx] += 1.0e-8;
    if(boundaries[mn] <= n && n <= boundaries[mx])
    {
        while(mx-mn > 1)
        {
            trial = (mn + mx)/2;
            if(n> boundaries[trial])
            {
                mn = trial;
            }
            else
            {
                mx = trial;
            }
        }
        return(mn);
    }
    else
    {
        return(-1);
    }
}

void GenerateCamRayIDPermutations(std::vector<std::vector<camrayid>> Lists, std::vector<std::vector<camrayid>>& result, unsigned int depth, std::vector<camrayid> current)
{
    if(depth == Lists.size())
    {
        result.push_back(current);
        return;
    }
    
    for(unsigned int i = 0; i < Lists[depth].size(); ++i)
    {
        std::vector<camrayid> tmp;
        tmp = current;
        tmp.push_back(Lists[depth][i]);
        GenerateCamRayIDPermutations(Lists, result, depth + 1, tmp);
    }
}

std::vector<traversedcell> DirectionalVoxelTraversal(transformedray ray, std::vector<std::vector<double>> bounds)
{
    std::vector<traversedcell> out;
    traversedcell tmptraversedcell;
    tmptraversedcell.camid = ray.camid;
    tmptraversedcell.rayid = ray.rayid;
    cellid curindex;
    curindex.xi = positionsorted(bounds[0], ray.x);
    curindex.yi = positionsorted(bounds[1], ray.y);
    curindex.zi = positionsorted(bounds[2], ray.z);
    
//    if(ray.camid == 1 && ray.rayid == 50)
//    {
//        std::cout << "\n raypos: = " << ray.x << " " << ray.y << " " << ray.z << "\n";
//        std::cout << "\ncurindex = " << curindex.xi << " " << curindex.yi << " " << curindex.zi << "\n";
//    }
    
    if(curindex.xi == -1 || curindex.yi == -1 || curindex.zi == -1)
    {
        std::cout << "ray starts outside bounds!\n";
        std::cout << "p:" << ray.x << ", " << ray.y <<  ", " << ray.z << "\n";
        std::cout << "v:" << ray.vx << ", " << ray.vy <<  ", " << ray.vz << "\n";
        std::cout << "v:" << curindex.xi << ", " << curindex.yi  <<  ", " << curindex.zi  << "\n";
        return(out);
    }
    else
    {
        long direction[3];
        direction[0] = Sgn(ray.vx);
        direction[1] = Sgn(ray.vy);
        direction[2] = Sgn(ray.vz);
        
        cellvisit tmpcellvisit;
        std::vector<cellvisit> cellvisits;
        for(int i = 0; i<3; i++)
        {
            double b = 0;
            double v = 0;
            switch (i) {
                case 0:
                    b = ray.x;
                    v = ray.vx;
                    break;
                case 1:
                    b = ray.y;
                    v = ray.vy;
                    break;
                case 2:
                    b = ray.z;
                    v = ray.vz;
            }
            for(unsigned int j = 0; j < bounds[i].size(); j++)
            {
                tmpcellvisit.dimension = i;
                tmpcellvisit.edgeQ = false;
                if(direction[i] == 1 && j == bounds[i].size() -1)
                {
                    tmpcellvisit.edgeQ = true;
                }
                if(direction[i] == -1 && j == 0)
                {
                    tmpcellvisit.edgeQ = true;
                }
                tmpcellvisit.time = specialdivision(bounds[i][j]-b,v);
                if(tmpcellvisit.time>1e-10)     // Only add those for which t is positive.  'ray is moving forward'
                {
                    cellvisits.push_back(tmpcellvisit);
                }
            }
        }
        
        std::sort(cellvisits.begin(), cellvisits.end(), comparecellvisitbytime);        // Sort by time of arrival
        
//       std::cout << "########## Edges " << ray.camid << "." << ray.rayid << " : ";
//        for(auto i: cellvisits)
//        {
//            std::cout << i.edgeQ << " " << i.time << " || ";
//        }
        
        for(long i = 0; cellvisits.size()-1;i++)
        {
            if(cellvisits[i].edgeQ == true)
            {
//                std::cout << "len = " << i << " \n";
                cellvisits.resize(i);
                break;
            }
        }
        
//        if(ray.camid == 1 && ray.rayid == 50)
//        {
//            std::cout << "Postlength: " << cellvisits.size() << " // ";
//            std::cout << "curindex " << curindex.xi << "."  << curindex.yi << "."  << curindex.zi << "\n";
//            for(auto cv:cellvisits)
//            {
//                std::cout << cv.dimension << " " <<cv.time << " " << cv.edgeQ << "\n";
//            }
//        }
        
        tmptraversedcell.cellid = curindex;     // Start at curindex
        out.push_back(tmptraversedcell);
        for (auto i: cellvisits)                // For each arrival time jump to next one...
        {
            switch (i.dimension) {
                case 0:
                    tmptraversedcell.cellid.xi += direction[i.dimension];
                    break;
                case 1:
                    tmptraversedcell.cellid.yi += direction[i.dimension];
                    break;
                case 2:
                    tmptraversedcell.cellid.zi += direction[i.dimension];
                    break;
                default:
                    break;
            }
            out.push_back(tmptraversedcell);
        }
        
//        if(ray.camid == 1 && ray.rayid == 50)
//        {
//            std::cout << "\n\n";
//            for(auto i: out)
//            {
//                std::cout << i.cellid.xi << " " << i.cellid.yi << " " << i.cellid.zi << "\n";
//            }
//            std::cout << "\n\n";
//        }
        
        // Use neighbours to 'expand out' out, and remove duplicates.
        long len = out.size();
        for(long i = 0; i <len; i++)
        {
            for(auto j: neighbours)
            {
                tmptraversedcell = out[i];
                tmptraversedcell.cellid.xi += j[0];
                tmptraversedcell.cellid.yi += j[1];
                tmptraversedcell.cellid.zi += j[2];
                out.push_back(tmptraversedcell);
            }
        }
        // Sort + delete duplicates
        std::sort(out.begin(),out.end(),comparetraversedcell);
        out.erase(unique(out.begin(), out.end(), comparetraversedcellsamecell), out.end());
        

        
        return(out);
    }
}

transformedray PrepareRay(ray r, boundingboxspec bb)
{
    transformedray out;
    out.camid = r.camid;
    out.rayid = r.rayid;
    out.x = r.x;
    out.y = r.y;
    out.z = r.z;
    double nrm = sqrt(r.vx*r.vx + r.vy*r.vy + r.vz*r.vz);
    out.vx = r.vx / nrm;
    out.vy = r.vy / nrm;
    out.vz = r.vz / nrm;


    if(bb.xmin < out.x && out.x < bb.xmax && bb.ymin < out.y && out.y < bb.ymax && bb.zmin < out.z && out.z < bb.zmax)
    {
        out.inside = true;
        out.hit = true;
        return out;
    }
    else
    {
        double t[6];
        t[0] = specialdivision(bb.xmin - out.x,out.vx);
        t[1] = specialdivision(bb.xmax - out.x,out.vx);
        t[2] = specialdivision(bb.ymin - out.y,out.vy);
        t[3] = specialdivision(bb.ymax - out.y,out.vy);
        t[4] = specialdivision(bb.zmin - out.z,out.vz);
        t[5] = specialdivision(bb.zmax - out.z,out.vz);
        
        std::vector<double> tmpip;
        std::vector<std::vector<double>> ip;
        double ti;
        for(int i=0; i<6; i++)
        {
            tmpip.clear();
            ti = t[i];
            if(ti == INFINITY || ti == -1*INFINITY)
            {
                tmpip.push_back(INFINITY);
                tmpip.push_back(INFINITY);
                tmpip.push_back(INFINITY);
            }
            else{
                tmpip.push_back(out.x + out.vx*ti);
                tmpip.push_back(out.y + out.vy*ti);
                tmpip.push_back(out.z + out.vz*ti);
            }
            ip.push_back(tmpip);
        }
        
        std::vector<bool> atfaces;
        atfaces.push_back(atface(bb.ymin, bb.ymax, ip[0][1]) && atface(bb.zmin, bb.zmax, ip[0][2]));
        atfaces.push_back(atface(bb.ymin, bb.ymax, ip[1][1]) && atface(bb.zmin, bb.zmax, ip[1][2]));
        atfaces.push_back(atface(bb.xmin, bb.xmax, ip[2][0]) && atface(bb.zmin, bb.zmax, ip[2][2]));
        atfaces.push_back(atface(bb.xmin, bb.xmax, ip[3][0]) && atface(bb.zmin, bb.zmax, ip[3][2]));
        atfaces.push_back(atface(bb.xmin, bb.xmax, ip[4][0]) && atface(bb.ymin, bb.ymax, ip[4][1]));
        atfaces.push_back(atface(bb.xmin, bb.xmax, ip[5][0]) && atface(bb.ymin, bb.ymax, ip[5][1]));
        std::vector<hitpoint> hitpoints;
        hitpoint tmphp;
        for(int i = 0; i<6; i++)
        {
            if(atfaces[i]){
                tmphp.t = t[i];
                tmphp.posx = ip[i][0];
                tmphp.posy = ip[i][1];
                tmphp.posz = ip[i][2];
                hitpoints.push_back(tmphp);
            }
        }
        
        if(hitpoints.size()>0)  // This should almost always be of length 2 if it hits the bb; unless it exactly hits the corner/edge of the bb...
        {
            std::sort(hitpoints.begin(), hitpoints.end(), comparehitpointbytime);   // Sort by arrival time
            out.hit = true;
            out.inside = false;
            out.x = hitpoints[0].posx;
            out.y = hitpoints[0].posy;
            out.z = hitpoints[0].posz;
            return(out);
        }
        else
        {
            out.hit = false;
            out.inside = false;
            return(out);
        }
    }
}

candidatematch ClosestPointToLines(std::map<std::pair<int, int>,transformedray>& raydb, std::vector<camrayid> crids)
{
    candidatematch out;
    out.camrayids = crids;
    long len = crids.size();
    CMatrix rhs(3,1);
    CMatrix lhs(3,3);
    for(int i=0;i<3;i++)            // Identity matrix and zero rhs / Fill in matrix for each 3 dimensions
    {
        rhs.m_pData[i][0] = 0.0;
        lhs.m_pData[i][i] = len;
    }
    
    for(int i=0;i<len;i++)            // Fill in rhs and lhs for each ray
    {
        transformedray r;
        r = raydb[std::make_pair(crids[i].camid,crids[i].rayid)];
        
        lhs.m_pData[0][0] -= r.vx*r.vx;
        lhs.m_pData[1][0] -= r.vx*r.vy;
        lhs.m_pData[2][0] -= r.vx*r.vz;
        lhs.m_pData[0][1] -= r.vy*r.vx;
        lhs.m_pData[1][1] -= r.vy*r.vy;
        lhs.m_pData[2][1] -= r.vy*r.vz;
        lhs.m_pData[0][2] -= r.vz*r.vx;
        lhs.m_pData[1][2] -= r.vz*r.vy;
        lhs.m_pData[2][2] -= r.vz*r.vz;
        
        double adotd = r.x*r.vx + r.y*r.vy + r.z*r.vz;
        rhs.m_pData[0][0] += r.x - r.vx*adotd;
        rhs.m_pData[1][0] += r.y - r.vy*adotd;
        rhs.m_pData[2][0] += r.z - r.vz*adotd;
    }

    // Solve it
    lhs = lhs.Inverse();
    CMatrix sol(3,1);
    sol = lhs*rhs;
    //std::cout << sol;
    double error = 0;
    double solminusax, solminusay, solminusaz;
    double tmpx,tmpy,tmpz;
    for(int i=0;i<len;i++)            // calculate distances
    {
        transformedray r;
        r = raydb[std::make_pair(crids[i].camid,crids[i].rayid)];
        
        solminusax = sol.m_pData[0][0] - r.x;
        solminusay = sol.m_pData[1][0] - r.y;
        solminusaz = sol.m_pData[2][0] - r.z;
        
        tmpx = r.vz * solminusay - r.vy * solminusaz;
        tmpy = solminusaz * r.vx - solminusax * r.vz;
        tmpz = solminusax * r.vy - solminusay * r.vx;
        
        error += tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
    }
    error = sqrt(error/len);
    
    out.matchx = sol.m_pData[0][0];
    out.matchy = sol.m_pData[1][0];
    out.matchz = sol.m_pData[2][0];
    out.matcherror = error;
    return(out);
}

std::vector<candidatematch> SpaceTraversalMatching(const std::vector<ray>& raydata, boundingboxspec bb, std::vector<std::vector<double>> bounds, int maxmatchesperray, unsigned int mincameras, double maxdistance, double multiplematchesperraymindistance)
{
    //std::cout << "Bounding box: " << bb.xmin << " - " << bb.xmax << " :: "  << bb.ymin << " - " << bb.ymax << " :: " << bb.zmin << " - " << bb.zmax << "\n";
    //std::cout << "First: " << raydata[0].camid << "-" << raydata[0].rayid << "\n";
    //std::cout << "Last: " << raydata[raydata.size()-1].camid << "-" << raydata[raydata.size()-1].rayid << "\n";
    
    // Prepare the rays
    std::map<std::pair<int, int>,transformedray> raydb;         // Store a dictionary of (cameraid, rayid): transformedray
    std::vector<transformedray> validrays;                      // Store the transformed rays
    std::map<int,int> numrays;                                  // Store the number of rays per camera
    std::map<int,int> invalidcounter;                           // Store the number that are invalid
    for(unsigned int i = 0; i < raydata.size(); i++)
    {
        if(numrays.find(raydata[i].camid) == numrays.end())     // Key does not yet exist  // maintain counter
            numrays.insert(std::make_pair(raydata[i].camid,1));
        else // Key exists
            numrays[raydata[i].camid]++;
        
        transformedray freshray;
        freshray = PrepareRay(raydata[i],bb);
        if(freshray.hit)
        {
            raydb.insert(std::make_pair(std::make_pair(freshray.camid, freshray.rayid),freshray));
            validrays.push_back(freshray);
        }
        else
        {
            if(invalidcounter.find(freshray.camid) == invalidcounter.end()) // Key does not yet exist  // maintain counter
                invalidcounter.insert(std::make_pair(freshray.camid,1));
            else // Key exists
                invalidcounter[freshray.camid]++;
        }
    }
    
    std::cout << "# of rays for each camera: {";
    for (auto kv = numrays.begin(); kv != numrays.end();) {
        std::cout << kv->first << ": " << kv->second;
        if(++kv != numrays.end())
        {
            std::cout << ", ";
        }
    }
    std::cout << "}\n";
    
    std::cout << "# of rays that miss the bounding box: {";
    for (auto kv = invalidcounter.begin(); kv != invalidcounter.end();) {
        std::cout << kv->first << ": " << kv->second;
        if(++kv != invalidcounter.end())
        {
            std::cout << ", ";
        }
    }
    std::cout << "}\n";
    
//    for (auto kv : invalidcounter) {
//        std::cout << "Camera " << kv.first << " has " << kv.second << " rays that miss the bounding box\n";
//    }
    
    // Do the traversing
    std::vector<traversedcell> traversed;
    std::vector<traversedcell> tmptraverse;
    std::vector<traversedcell> tmptraverse2;
    transformedray revvr;
    //validrays.resize(48);
    for (transformedray vr : validrays)
    {
        tmptraverse = DirectionalVoxelTraversal(vr,bounds);         // Traverse forwards
        if(vr.inside)                                               // Ray is inside, traverse also backwards
        {
            revvr = vr;
            revvr.vx *= -1;
            revvr.vy *= -1;
            revvr.vz *= -1;
            tmptraverse2 = DirectionalVoxelTraversal(revvr,bounds); // Combine these two…
            tmptraverse.insert(tmptraverse.end(), tmptraverse2.begin(), tmptraverse2.end()); //
        }
        std::sort(tmptraverse.begin(),tmptraverse.end(),comparetraversedcell);
        tmptraverse.erase(std::unique(tmptraverse.begin(),tmptraverse.end(),comparetraversedcellsamecell),tmptraverse.end());
        //std::cout << tmptraverse.size() << " ";
        traversed.insert(traversed.end(), tmptraverse.begin(), tmptraverse.end());
    }
//    for(auto i: traversed)
//    {
//            std::cout << i.camid << "." << i.rayid << "  " << i.cellid.xi << "," << i.cellid.yi << "," << i.cellid.xi << "\n";
//    }
    
    std::cout << "# of voxels traversed after expansion: " << traversed.size() << "\n";
    std::sort(traversed.begin(),traversed.end(),comparetraversedcell);
    
    std::vector<groupedcell> groupedcells;
    groupedcell tmpgroupedcell;
    tmpgroupedcell.camrayids.push_back({traversed[0].camid, traversed[0].rayid});
    tmpgroupedcell.cellid = traversed[0].cellid;
    for(unsigned long i = 1; i < traversed.size(); i++)
    {

        if(tmpgroupedcell.cellid.xi == traversed[i].cellid.xi && tmpgroupedcell.cellid.yi == traversed[i].cellid.yi && tmpgroupedcell.cellid.zi == traversed[i].cellid.zi)
        {
            tmpgroupedcell.camrayids.push_back({traversed[i].camid, traversed[i].rayid});
        }
        else
        {
            if(tmpgroupedcell.camrayids.size() >= mincameras)       // Rough prune
            {
                std::vector<camrayid> tmpcamrayids = tmpgroupedcell.camrayids;
                std::sort(tmpcamrayids.begin(), tmpcamrayids.end(), comparecamrayidscam);
                long uniquecount = std::unique(tmpcamrayids.begin(), tmpcamrayids.end(), comparecamrayidsamecam) - tmpcamrayids.begin();
                if(uniquecount >= mincameras)                       // More careful prune
                {
                    groupedcells.push_back(tmpgroupedcell);
                }
            }
            tmpgroupedcell.cellid = traversed[i].cellid;
            tmpgroupedcell.camrayids.clear();
            tmpgroupedcell.camrayids.push_back({traversed[i].camid, traversed[i].rayid});
        }
    }
    // Process remainder:
    if(tmpgroupedcell.camrayids.size() >= mincameras)       // Rough prune
    {
        std::vector<camrayid> tmpcamrayids = tmpgroupedcell.camrayids;
        std::sort(tmpcamrayids.begin(), tmpcamrayids.end(), comparecamrayidscam);
        long uniquecount = std::unique(tmpcamrayids.begin(), tmpcamrayids.end(), comparecamrayidsamecam) - tmpcamrayids.begin();
        if(uniquecount >= mincameras)                       // More careful prune
        {
            groupedcells.push_back(tmpgroupedcell);
        }
    }
    std::cout << "Prune based on number of cameras: " << groupedcells.size() << "\n";
    
    // Get rid of cellID
    std::vector<std::vector<camrayid>> candidatepairs;
    std::vector<camrayid> tmpcamrayids;
    for(auto i: groupedcells)
    {
        tmpcamrayids = i.camrayids;
        std::sort(tmpcamrayids.begin(),tmpcamrayids.end(),comparecamrayidscamray);
        candidatepairs.push_back(tmpcamrayids);
        
    }
    // This most of the time semi-sorted, so we can unique first and remove a lot…
    candidatepairs.erase(std::unique(candidatepairs.begin(), candidatepairs.end(), comparecamrayidsidentical), candidatepairs.end());
    std::sort(candidatepairs.begin(),candidatepairs.end(),comparecamrayidsordered);         // Now sort followed by another unique…
    candidatepairs.erase(std::unique(candidatepairs.begin(), candidatepairs.end(), comparecamrayidsidentical), candidatepairs.end());
    
    std::vector<std::vector<camrayid>> candidates;
    for(std::vector<camrayid> cand: candidatepairs)
    {
        // Split candidates in to groups based on camera:
        std::vector<std::vector<camrayid>> groupedcandidates;
        std::vector<camrayid> tmpgroupedcandidates;
        tmpgroupedcandidates.push_back(cand[0]);
        for(unsigned long i = 1; i < cand.size(); i++)
        {
            
            if(tmpgroupedcandidates[0].camid == cand[i].camid)
            {
                tmpgroupedcandidates.push_back(cand[i]);
            }
            else
            {
                groupedcandidates.push_back(tmpgroupedcandidates);
                tmpgroupedcandidates.clear();
                tmpgroupedcandidates.push_back(cand[i]);
            }
        }
        groupedcandidates.push_back(tmpgroupedcandidates);
        std::vector<camrayid> tmp;
        std::vector<std::vector<camrayid>> newcandidates;
        
        GenerateCamRayIDPermutations(groupedcandidates, newcandidates, 0, tmp);
        candidates.insert(candidates.end(), newcandidates.begin(), newcandidates.end());
    }
//    std::cout << "Candidates size: " << candidates.size() << "\n";

    std::sort(candidates.begin(),candidates.end(),comparecamrayidsordered);
    candidates.erase(std::unique(candidates.begin(), candidates.end(), comparecamrayidsidentical), candidates.end());
    
    std::cout << "Duplicate candidates removed: " << candidates.size() << "\n";
//    for(auto i: candidates)
//    {
//        std::cout << "candidates: ";
//        for (auto j: i)
//        {
//            std::cout << j.camid << "." <<j .rayid << " ";
//        }
//        std::cout << "\n";
//    }
    
    // Calculate match positions and errors for each of the candidates
    std::vector<candidatematch> candidatematches;
    for(auto cand: candidates)
    {
        candidatematch tmp;
        tmp = ClosestPointToLines(raydb,cand);
        candidatematches.push_back(tmp);
    }
    
    std::sort(candidatematches.begin(), candidatematches.end(), comparecandidatematches);
    
    
//    for(auto i: candidatematches)
//    {
//        std::cout << "Candidates match: ";
//        for (auto j: i.camrayids)
//        {
//            std::cout << j.camid << "." <<j .rayid << "\t";
//        }
//        std::cout << "\t| "  << i.matcherror << "\t| " << i.matchx << "\t"  << i.matchy << "\t"  << i.matchz << "\t" <<  "\n";
////        for (auto j: i.camrayids)
////        {
////            transformedray r;
////            r = raydb[std::make_pair(j.camid,j.rayid)];
////            std::cout << r.x << "\t" << r.y << "\t" << r.z << "\t"  << r.vx << "\t" << r.vy << "\t" << r.vz << "\t\n";
////        }
////        std::cout << "\n";
//    }
    
    // Select best matches from candidates
    std::cout << "Selecting the best matches with upto " << maxmatchesperray << " match(es)/ray out of " << candidates.size() <<" candidates\n";
    std::vector<candidatematch> approvedmatches;
    std::map<std::pair<long, long>,int> matchcounter;
    std::map<std::pair<long, long>,std::vector<candidatematch>> approvedmatches_per_ray;
    unsigned int removed_because_sphere = 0;
    approvedmatches.clear();
    matchcounter.clear();
    bool valid;
    int matches;
    for(candidatematch cand: candidatematches)
    {
        if(cand.matcherror < maxdistance)
        {
            valid = true;
            for(auto camrayid: cand.camrayids)
            {
                std::pair<long, long> idpair;
                idpair = std::make_pair(camrayid.camid,camrayid.rayid);
                if(!(matchcounter.find(idpair) == matchcounter.end())) // If it exists
                {
                    // If there are more than maxmatchesperray for this ray, valid=false
                    matches = matchcounter[idpair];
                    if(matches >= maxmatchesperray)
                    {
                        valid = false;
                        break;
                    }

                    // If there is another match for this ray in the exclusion sphere, valid=false
                    for(auto othermatch: approvedmatches_per_ray[idpair])
                    {
                        double distance = sqrt(pow(cand.matchx-othermatch.matchx, 2) + pow(cand.matchy-othermatch.matchy, 2) + pow(cand.matchz-othermatch.matchz, 2));
                        if (distance < multiplematchesperraymindistance)
                        {
                            valid = false;
                            removed_because_sphere++;
                            break;
                        }
                    }
                    if (!valid)
                    {
                        break;
                    }
                    
                }
            }
            if(valid)
            {
                // Add to approved matches
                for(auto camrayid: cand.camrayids)  // Add match counter
                {
                    std::pair<long, long> idpair;
                    idpair = std::make_pair(camrayid.camid,camrayid.rayid);
                    if(matchcounter.find(idpair) == matchcounter.end())                     // If does not exists
                    {
                       matchcounter.insert(std::pair<std::pair<long,long>,int>(idpair,1));  // Add to counter
                       std::vector<candidatematch> v;
                       v.push_back(cand);
                       approvedmatches_per_ray.insert(std::pair<std::pair<long,long>, std::vector<candidatematch>>(idpair, v));
                    }
                    else
                    {
                        matchcounter[idpair]++;                                            // Increment counter
                        approvedmatches_per_ray[idpair].push_back(cand);
                    }
                }
                approvedmatches.push_back(cand);
            }
        }
    }
    
    std::cout << "Selecting done. " << approvedmatches.size() << " matched found (out of " << candidates.size() << " candidates)\n";
    if (multiplematchesperraymindistance > 0.0)
    {
        std::cout << "Matches remove because of minimum distance for multiple matches per ray " << removed_because_sphere << std::endl;
        std::cout << "Minimum distance for multiple matches per ray is " << multiplematchesperraymindistance << std::endl;
    }
    else
    {
        if (removed_because_sphere != 0)
        {
            std::cout << "****** BUG ******: removed_because_sphere should be 0 !!" << std::endl;
        }
    }
    
//    for(auto i: approvedmatches)
//    {
//        std::cout << "Approved match: ";
//        for (auto j: i.camrayids)
//        {
//            std::cout << j.camid << "." <<j .rayid << "\t";
//        }
//        std::cout << "\t| "  << i.matcherror << "\n";
//    }
    
    return(approvedmatches);
}

void init(){
    neighbours.clear();                 // Presumably it is already empty…
    neighbours.push_back({-1,0,0});     // Add neighbours that are above/below/left/right/front/back
    neighbours.push_back({0,-1,0});
    neighbours.push_back({0,0,-1});
    neighbours.push_back({0,0,1});
    neighbours.push_back({0,1,0});
    neighbours.push_back({1,0,0});
}
