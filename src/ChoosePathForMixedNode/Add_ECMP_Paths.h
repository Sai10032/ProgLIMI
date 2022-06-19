#ifndef ADD_ECMP_PATHS_H
#define ADD_ECMP_PATHS_H
#include <math.h>
#include <algorithm>
#include <fstream>
#include <map>
#include "Boost_Graph.h"
void Compute_ECMP_Paths_From_Between_Two_Nodes(BGraph& bgraph, unsigned s, unsigned d, list<PathType>& paths);

#endif
