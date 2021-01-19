#include <vector>
#include <set>

#include "weights/VoronoiUtils.h"
#include "weights/PolysToContigWeights.h"
#include "weights/GalWeight.h"
#include "weights/GeodaWeight.h"
#include "SpatialIndAlgs.h"
#include "gda_interface.h"
#include "gda_weights.h"

GeoDaWeight* contiguity_weights(bool is_queen,
                                AbstractGeoDa* geoda,
                                unsigned int order,
                                bool include_lower_order,
                                double precision_threshold)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();
    GalWeight* poW = new GalWeight;
    poW->num_obs = num_obs;
    poW->is_symmetric = true;
    poW->symmetry_checked = true;

    if (geoda->GetMapType() == gda::POINT_TYP) {
        std::vector<std::set<int> > nbr_map;
        const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();
        std::vector<double> x(num_obs), y(num_obs);
        for (int i=0; i<num_obs; ++i) {
            x[i] = centroids[i]->x;
            y[i] = centroids[i]->y;
        }
        Gda::VoronoiUtils::PointsToContiguity(x, y, is_queen, nbr_map);
        poW->gal = Gda::VoronoiUtils::NeighborMapToGal(nbr_map);
        //gda::PointsToContiguity(x, y, is_queen, nbr_map);
        //poW->gal = Gda::NeighborMapToGal(nbr_map);
        if (order > 1) {
            Gda::MakeHigherOrdContiguity(order, num_obs, poW->gal, include_lower_order);
        }

    } else if (geoda->GetMapType() == gda::POLYGON) {
        poW->gal = PolysToContigWeights(geoda->GetMainMap(), is_queen, precision_threshold);
        if (order > 1) {
            Gda::MakeHigherOrdContiguity(order, num_obs, poW->gal, include_lower_order);
        }

    } else {
        // line_type not supported yet, should be detected at script side
        delete poW;
        return 0;
    }

    poW->GetNbrStats();
    return (GeoDaWeight*)poW;
}

GeoDaWeight* gda_queen_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = true;
    return contiguity_weights(is_queen, geoda, order, include_lower_order, precision_threshold);
}

GeoDaWeight* gda_rook_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = false;
    return contiguity_weights(is_queen, geoda, order, include_lower_order, precision_threshold);
}

GeoDaWeight* gda_knn_weights(AbstractGeoDa* geoda, unsigned int k,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagnals,
                             const std::string& polyid)
{
    if (geoda == 0) return 0;

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    GwtWeight* poW = SpatialIndAlgs::knn_build(centroids, k, is_arc,
                                               is_mile, is_inverse, power,
                                               kernel, bandwidth, adaptive_bandwidth, use_kernel_diagnals);
    poW->GetNbrStats();
    return (GeoDaWeight*)poW;
}

double gda_min_distthreshold(AbstractGeoDa* geoda, bool is_arc, bool is_mile)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    std::vector<double> x(num_obs), y(num_obs);
    for (int i=0; i<num_obs; ++i) {
        x[i] = centroids[i]->x;
        y[i] = centroids[i]->y;
    }

    double max_1nn_dist = SpatialIndAlgs::find_max_1nn_dist(x, y, is_arc, is_mile);
    return max_1nn_dist;
}

GeoDaWeight* gda_distance_weights(AbstractGeoDa* geoda, double dist_thres,
                                  const std::string& polyid,
                                  double power,
                                  bool is_inverse,
                                  bool is_arc,
                                  bool is_mile,
                                  const std::string& kernel,
                                  bool use_kernel_diagnals)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();
    
    std::vector<double> x(num_obs), y(num_obs);
    for (int i=0; i<num_obs; ++i) {
        x[i] = centroids[i]->x;
        y[i] = centroids[i]->y;
    }
    dist_thres = dist_thres * 1.00000; //m_thres_delta_factor
    GwtWeight* poW = SpatialIndAlgs::thresh_build(x, y, dist_thres, power, is_arc, is_mile,
                                                  kernel, use_kernel_diagnals);

    poW->GetNbrStats();
    return (GeoDaWeight*)poW;
}
