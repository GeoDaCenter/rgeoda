//
// Created by Xun Li on 1/14/21.
//

#ifndef GEODA_SCHC_WRAPPER_H
#define GEODA_SCHC_WRAPPER_H

#include <vector>
#include <string>

class GeoDaWeight;

class schc_wrapper {
public:
    schc_wrapper(unsigned int k,
                   GeoDaWeight *w,
                   const std::vector<std::vector<double> > &data,
                   unsigned int redcap_method,
                   const std::string &distance_method,
                   const std::vector<double>& bound_vals,
                   double min_bound);

    virtual ~schc_wrapper();

    const std::vector<int> GetFlatClusters();

    const std::vector<std::vector<int> > GetClusters();

private:
    int num_obs;

    std::vector<std::vector<int> > cluster_ids;
};


#endif //GEODA_SCHC_WRAPPER_H
