#ifndef __JSGEODSA_GDA_WEIGHTS__
#define __JSGEODSA_GDA_WEIGHTS__

#include <string>

class AbstractGeoDa;

class GeoDaWeight;

// APIs of weights creation
/**
 *
 * @param geoda
 * @param polyid
 * @param order
 * @param include_lower_order
 * @param precision_threshold
 * @return
 */
GeoDaWeight* gda_queen_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold);

GeoDaWeight* gda_rook_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold);

GeoDaWeight* gda_knn_weights(AbstractGeoDa* geoda, unsigned int k,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagnals,
                             const std::string& polyid);

double gda_min_distthreshold(AbstractGeoDa* geoda, bool is_arc, bool is_mile);

GeoDaWeight* gda_distance_weights(AbstractGeoDa* geoda, double dist_thres,
                                  const std::string& polyid,
                                  double power,
                                  bool is_inverse,
                                  bool is_arc,
                                  bool is_mile,
                                  const std::string& kernel,
                                  bool use_kernel_diagonals);
#endif
