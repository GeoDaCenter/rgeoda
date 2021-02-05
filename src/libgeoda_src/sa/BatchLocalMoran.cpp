//
// Created by Xun Li on 2019-06-05.
//

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"
#include "BatchLocalMoran.h"


BatchLocalMoran::~BatchLocalMoran() {

}

BatchLocalMoran::BatchLocalMoran(int num_obs,
                 GeoDaWeight *w,
                 const std::vector<std::vector<double> > &_data,
                 const std::vector<std::vector<bool> > &_undefs,
                 double _significance_cutoff,
                 int _nCPUs, int _perm, uint64_t _last_seed)
: BatchLISA(num_obs, w, _undefs, _significance_cutoff, _nCPUs, _perm, _last_seed),
CLUSTER_NOT_SIG(0),
CLUSTER_HIGHHIGH(1),
CLUSTER_LOWLOW(2),
CLUSTER_LOWHIGH(3),
CLUSTER_HIGHLOW(4),
CLUSTER_UNDEFINED(5),
CLUSTER_NEIGHBORLESS(6),
data(_data)
{
    labels.push_back("Not significant");
    labels.push_back("High-High");
    labels.push_back("Low-Low");
    labels.push_back("High-Low");
    labels.push_back("Low-High");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#FF0000");
    colors.push_back("#0000FF");
    colors.push_back("#a7adf9");
    colors.push_back("#f4ada8");
    colors.push_back("#464646");
    colors.push_back("#999999");

    num_batch = data.size();
    for (int i=0; i<num_batch; ++i) {
        GenUtils::StandardizeData(data[i], undefs[i]);
    }

    Run();
}

void BatchLocalMoran::ComputeLoalSA() {
    for (int v=0; v<num_batch; ++v) {
        for (int i = 0; i < num_obs; i++) {
            if (undefs[v][i]) {
                lag_vec[v][i] = 0;
                lisa_vec[v][i] = 0;
                cluster_vec[v][i] = CLUSTER_UNDEFINED;

            } else {
                if (weights->GetNbrSize(i) == 0) {
                    cluster_vec[v][i] = CLUSTER_NEIGHBORLESS;
                } else {
                    double sp_lag = 0;
                    const std::vector<long> &nbrs = weights->GetNeighbors(i);
                    unsigned int nn = 0;
                    for (size_t j = 0; j < nbrs.size(); ++j) {
                        if (nbrs[j] != i && !undefs[v][nbrs[j]]) { // not including the value at the location
                            sp_lag += data[v][nbrs[j]];
                            nn += 1;
                        }
                    }
                    sp_lag = sp_lag / nn;
                    lag_vec[v][i] = sp_lag;
                    lisa_vec[v][i] = data[v][i] * sp_lag;
                    // assign the cluster
                    if (data[v][i] > 0 && sp_lag < 0) cluster_vec[v][i] = CLUSTER_HIGHLOW;
                    else if (data[v][i] < 0 && sp_lag > 0) cluster_vec[v][i] = CLUSTER_LOWHIGH;
                    else if (data[v][i] < 0 && sp_lag < 0) cluster_vec[v][i] = CLUSTER_LOWLOW;
                    else cluster_vec[v][i] = CLUSTER_HIGHHIGH; //data1[i] > 0 && Wdata > 0

                }
            }
        }
    }
}

void BatchLocalMoran::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors,
        std::vector<std::vector<double> >& permutedSA)
{
    for (int v=0; v<num_batch; ++v) {
        int validNeighbors = 0;
        double permutedLag = 0;
        int numNeighbors = permNeighbors.size();
        // use permutation to compute the lag
        // compute the lag for binary weights
        for (int cp = 0; cp < numNeighbors; cp++) {
            int nb = permNeighbors[cp];
            if (!undefs[v][nb]) {
                permutedLag += data[v][nb];
                validNeighbors++;
            }
        }
        //NOTE: we shouldn't have to row-standardize or
        // multiply by data1[cnt]
        if (validNeighbors > 0 && row_standardize) {
            permutedLag /= validNeighbors;
        }
        const double localMoranPermuted = permutedLag * data[v][cnt];
        permutedSA[v][perm] = localMoranPermuted;
    }
}

std::vector<uint64_t> BatchLocalMoran::CountLargerSA(int cnt, const std::vector<std::vector<double> >& permutedSA)
{
    std::vector<uint64_t> results(num_batch);
    for (int v=0; v<num_batch; ++v) {
        uint64_t countLarger = 0;
        for (int i = 0; i < permutations; ++i) {
            if (permutedSA[v][i] >= lisa_vec[v][cnt]) {
                countLarger += 1;
            }
        }

        // pick the smallest counts
        if (permutations - countLarger <= countLarger) {
            countLarger = permutations - countLarger;
        }
        results[v] = countLarger;
    }
    return results;
}

std::vector<int> BatchLocalMoran::GetClusterIndicators(int idx)
{
    std::vector<int> clusters(num_obs);

    double cuttoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if (sig_local_vec[idx][i] > cuttoff &&
                (const unsigned long)cluster_vec[idx][i] != CLUSTER_UNDEFINED &&
                (const unsigned long)cluster_vec[idx][i] != CLUSTER_NEIGHBORLESS)
        {
            clusters[i] = CLUSTER_NOT_SIG;
        } else {
            clusters[i] = cluster_vec[idx][i];
        }
    }
    return clusters;
}

