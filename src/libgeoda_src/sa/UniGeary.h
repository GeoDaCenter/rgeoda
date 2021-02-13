//
// Created by Xun Li on 9/27/19.
//

#ifndef GEODA_UNIGEARY_H
#define GEODA_UNIGEARY_H

#include <vector>
#include "LISA.h"

class UniGeary : public LISA {
    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_HIGHHIGH;
    const unsigned long CLUSTER_LOWLOW;
    const unsigned long CLUSTER_OTHERPOS;
    const unsigned long CLUSTER_NEGATIVE;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    UniGeary(int num_obs,
             GeoDaWeight* w,
             const std::vector<double>& data,
             const std::vector<bool>& undefs,
             double significance_cutoff,
             int nCPUs, int permutations,
             const std::string& _permutation_method,
             uint64_t last_seed_used);


    virtual ~UniGeary();

    virtual void ComputeLoalSA() ;

    virtual void PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA);

    virtual void PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors, std::vector<double>& permutedSA){}

    virtual uint64_t CountLargerSA(int cnt, const std::vector<double>& permutedSA);

    virtual std::vector<int> GetClusterIndicators();

protected:
    std::vector<double> data;
    std::vector<double> data_square;
};



#endif //GEODA_UNIGEARY_H
