/*
 * Main STM functions
 *
 */

#include <vector>
#include "STM_types.h"

bool comparecellvisitbytime(const cellvisit &a, const cellvisit &b);
bool comparecandidatematches(const candidatematch &a, const candidatematch &b);
bool comparecamrayidscam(const camrayid &a, const camrayid &b);
bool comparecamrayididentical(const camrayid &a, const camrayid &b);
bool comparecamrayidscamray(const camrayid &a, const camrayid &b);
bool comparecamrayidsamecam(const camrayid &a, const camrayid &b);
bool comparecamrayidsidentical(const std::vector<camrayid> &a, const std::vector<camrayid> &b);
bool comparecamrayidsordered(const std::vector<camrayid> &a, const std::vector<camrayid> &b);
bool comparetraversedcell(const traversedcell &a, const traversedcell &b);
bool comparetraversedcellsamecell(const traversedcell &a, const traversedcell &b);
bool comparehitpointbytime(const hitpoint &a, const hitpoint &b);
long positionsorted(std::vector<double> boundaries,double n);

void GenerateCamRayIDPermutations(std::vector<std::vector<camrayid>> Lists, std::vector<std::vector<camrayid>>& result, unsigned int depth, std::vector<camrayid> current);
std::vector<traversedcell> DirectionalVoxelTraversal(transformedray ray, std::vector<std::vector<double>> bounds);
transformedray PrepareRay(ray r, boundingboxspec bb);
candidatematch ClosestPointToLines(std::map<std::pair<int, int>,transformedray>& raydb, std::vector<camrayid> crids);
std::vector<candidatematch> SpaceTraversalMatching(const std::vector<ray>& raydata, boundingboxspec bb, std::vector<std::vector<double>> bounds, int maxmatchesperray, unsigned int mincameras, double maxdistance, double multiplematchesperraymindistance);

void init();
