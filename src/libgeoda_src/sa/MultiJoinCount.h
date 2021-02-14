//
// Created by Xun Li on 2019-12-01.
//

#ifndef GEODA_MULTIJOINCOUNT_H
#define GEODA_MULTIJOINCOUNT_H

#include "LISA.h"

class MultiJoinCount : public LISA {
    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_SIG;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    MultiJoinCount(int num_obs, GeoDaWeight* w,
                   const std::vector<std::vector<double> >& data,
                   const std::vector<std::vector<bool> >& undefs,
                   double significance_cutoff,
                   int nCPUs, int permutations,
                   const std::string& _permutation_method,
                   uint64_t last_seed_used);

    virtual ~MultiJoinCount();

    virtual void CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start);

    virtual void PermCalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start);

    virtual void ComputeLoalSA() ;

    virtual void PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA);

    virtual void PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors, std::vector<double>&
            permutedSA);

    virtual uint64_t CountLargerSA(int cnt, const std::vector<double>& permutedSA);

    virtual std::vector<int> GetClusterIndicators();

protected:
    int num_vars;
    std::vector<int> zz;
    std::vector<std::vector<double> > data;
};

#endif //GEODA_MULTIJOINCOUNT_H
