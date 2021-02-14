//
// Created by Xun Li on 10/2/19.
//

#ifndef GEODA_UNIGSTAR_H
#define GEODA_UNIGSTAR_H


#include <vector>

#include "LISA.h"

class GeoDaWeight;

class UniGstar: public LISA {
    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_HIGHHIGH;
    const unsigned long CLUSTER_LOWLOW;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    UniGstar(int num_obs,
             GeoDaWeight* w,
             const std::vector<double>& data,
             const std::vector<bool>& undefs,
             double significance_cutoff,
             int nCPUs, int permutations,
             const std::string& _permutation_method,
             uint64_t last_seed_used);

    virtual ~UniGstar();

    virtual void ComputeLoalSA() ;

    virtual void PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA);

    virtual void PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors,
            std::vector<double>& permutedSA);

    virtual uint64_t CountLargerSA(int cnt, const std::vector<double>& permutedSA);

    virtual std::vector<int> GetClusterIndicators();

protected:
    std::vector<double> data;

    double sum_x;
    std::vector<bool> Gstar_defined;
};


#endif //GEODA_UNIGSTAR_H
