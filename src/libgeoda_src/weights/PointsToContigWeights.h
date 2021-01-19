#ifndef __JSGEODA_POINTSTOCONTIGWEIGHTS__
#define __JSGEODA_POINTSTOCONTIGWEIGHTS__

#include <vector>
#include <set>

namespace gda {
    bool PointsToContiguity(const std::vector<double>& x,
							const std::vector<double>& y,
							bool queen, // if false, then rook only
							std::vector<std::set<int> >& nbr_map);
}

#endif
