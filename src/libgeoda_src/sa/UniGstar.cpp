//
// Created by Xun Li on 10/2/19.
//
#include <math.h>
#include <vector>

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"

#include "UniGstar.h"

UniGstar::~UniGstar() {

}

UniGstar::UniGstar(int num_obs,
           GeoDaWeight *w,
           const std::vector<double> &_data,
           const std::vector<bool> &_undefs,
                   int _nCPUs, int _perm, uint64_t _last_seed)
        : LISA(num_obs, w, _undefs, _nCPUs, _perm, _last_seed), 
          CLUSTER_NOT_SIG(0),
          CLUSTER_HIGHHIGH(1),
          CLUSTER_LOWLOW(2),
          CLUSTER_UNDEFINED(3),
          CLUSTER_NEIGHBORLESS(4),
          data(_data), undefs(_undefs),
          sum_x(0)

{
    labels.push_back("Not significant");
    labels.push_back("High-High");
    labels.push_back("Low-Low");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#FF0000");
    colors.push_back("#0000FF");
    colors.push_back("#464646");
    colors.push_back("#999999");

    Gstar_defined.resize(num_obs, true);

    for (int i=0; i<num_obs; i++) {
        if (!undefs[i])
            sum_x += data[i];
    }

    Run();
}

void UniGstar::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i]) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;

        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                double lag = 0;
                const std::vector<long>& nbrs = weights->GetNeighbors(i);
                unsigned int nn = 0;
                for (size_t j=0; j<nbrs.size(); ++j) {
                    if (nbrs[j] != i &&  !undefs[nbrs[j]]) { // not including the value at the location
                        lag += data[ nbrs[j] ];
                        nn += 1;
                    }
                }
                // including self for G*
                lag += data[i];
                nn += 1;
                // row-standardize
                lag = lag / nn;
                lisa_vec[i] = lag / sum_x;
            }
        }
    }
    // mean G value
    unsigned int ng = 0;
    double mean_g = 0;
    for (int i=0; i<num_obs; ++i) {
        if (weights->GetNbrSize(i) == 0 || undefs[i] || Gstar_defined[i] == false)
            continue;
        mean_g += lisa_vec[i];
        ng += 1;
    }
    mean_g = mean_g / ng;

    // assign cluster
    for (int i=0; i<num_obs; ++i) {
        if (weights->GetNbrSize(i) == 0 || undefs[i] || Gstar_defined[i] == false)
            continue;

        if (lisa_vec[i] >= mean_g) {
            cluster_vec[i] = CLUSTER_HIGHHIGH;
        } else {
            cluster_vec[i] = CLUSTER_LOWLOW;
        }
    }
}

void UniGstar::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors,
                       std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    int numNeighbors = permNeighbors.size();

    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb] && nb != cnt) {
            permutedLag += data[nb];
            validNeighbors ++;
        }
    }
    // including self
    permutedLag += data[cnt];
    validNeighbors += 1;

    double permutedG = 0;
    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors;
        permutedG = permutedLag / sum_x;
    }
    permutedSA[perm] = permutedG;
}

uint64_t UniGstar::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
{
    uint64_t countLarger = 0;
    for (int i=0; i<permutations; ++i) {
        if (permutedSA[i] >= lisa_vec[cnt]) {
            countLarger += 1;
        }
    }

    // pick the smallest counts
    if (permutations-countLarger < countLarger) {
        countLarger = permutations-countLarger;
    }
    return countLarger;
}

std::vector<int> UniGstar::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cutoff = GetSignificanceCutoff();

    for (int i=0; i<num_obs; i++) {
        if ((const unsigned long)cluster_vec[i] == CLUSTER_UNDEFINED &&
                (const unsigned long)cluster_vec[i] == CLUSTER_NEIGHBORLESS)
            continue;

        if (sig_local_vec[i] > cutoff) {
            clusters[i] = CLUSTER_NOT_SIG;
        } else {
            clusters[i] = cluster_vec[i];
        }
    }
    return clusters;
}

