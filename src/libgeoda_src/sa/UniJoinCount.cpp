//
// Created by Xun Li on 9/27/19.
//
#include <vector>

#include "UniJoinCount.h"
#include "../weights/GeodaWeight.h"
#include "../GenUtils.h"

UniJoinCount::UniJoinCount(int num_obs, GeoDaWeight *w,
        const std::vector<double> &_data,
        const std::vector<bool> &_undefs,
        int _nCPUs, int _perm, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, _nCPUs, _perm, _last_seed), 
  CLUSTER_NOT_SIG(0),
  CLUSTER_SIG(1),
  CLUSTER_UNDEFINED(2),
  CLUSTER_NEIGHBORLESS(3),
  data(_data), undefs(_undefs)
{
    labels.push_back("Not significant");
    labels.push_back("Significant");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#348124"); // green
    colors.push_back("#464646");
    colors.push_back("#999999");

    Run();
}

UniJoinCount::~UniJoinCount() {

}

void UniJoinCount::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i]) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;
        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                if (data[i] > 0) { // x_i = 1
                    int nbr_size = weights->GetNbrSize(i);
                    const std::vector<long>& nbrs = weights->GetNeighbors(i);
                    for (int j=0; j<nbr_size; ++j) {
                        if (nbrs[j] != i &&  !undefs[nbrs[j]]) { // not including the value at the location
                            lisa_vec[i] += data[nbrs[j]];
                        }
                    }
                }
            }
        }
    }
}

void
UniJoinCount::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double> &permutedSA) {

    int validNeighbors = 0;
    double permutedLag = 0;
    int numNeighbors = permNeighbors.size();
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb]) {
            permutedLag += data[nb];
            validNeighbors ++;
        }
    }
    permutedSA[perm] = permutedLag;
}

uint64_t UniJoinCount::CountLargerSA(int cnt, const std::vector<double> &permutedSA) {
    uint64_t countLarger = 0;
    for (int i=0; i<permutations; ++i) {
        if (permutedSA[i] >= lisa_vec[cnt]) {
            countLarger += 1;
        }
    }

    // pick the smallest counts
    if (permutations-countLarger <= countLarger) {
        countLarger = permutations-countLarger;
    }
    return countLarger;
}

std::vector<int> UniJoinCount::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cuttoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if (sig_local_vec[i] <= cuttoff ) {
            if (lisa_vec[i] == 0) {
                clusters[i] = CLUSTER_NOT_SIG;
            } else {
                clusters[i] = CLUSTER_SIG;
            }
        } else {
            clusters[i] = CLUSTER_NOT_SIG;
        }
    }
    return clusters;
}

