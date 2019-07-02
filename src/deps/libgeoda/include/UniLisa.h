//
// Created by Xun Li on 2019-06-05.
//

#ifndef GEODA_UNILISA_H
#define GEODA_UNILISA_H

#include <vector>

#include "AbstractLocalSA.h"

class GeoDaWeight;

class UniLisa : public AbstractLocalSA {

    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_HIGHHIGH;
    const unsigned long CLUSTER_LOWLOW;
    const unsigned long CLUSTER_LOWHIGH;
    const unsigned long CLUSTER_HIGHLOW;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    UniLisa(int num_obs,
            const std::vector<double>& data,
            const std::vector<bool>& undefs,
            GeoDaWeight* w);


    virtual ~UniLisa();

    virtual void ComputeLoalSA() ;

    virtual bool LargerPermLocalSA(int cnt, std::vector<int> &permNeighbors);

    virtual std::vector<int> GetClusterIndicators();

    const std::vector<double> &GetLagValues() const;

    const std::vector<double> &GetLocalMoranValues() const;

protected:
    std::vector<double> data;
    std::vector<bool> undefs;

    std::vector<double> lag_vec;
    std::vector<double> localMoran_vec;


};

#endif //GEODA_UNILISA_H
