//
// Created by Xun Li on 9/27/19.
//

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"
#include "UniGeary.h"

UniGeary::UniGeary(int num_obs, GeoDaWeight *w, const std::vector<double> &_data, const std::vector<bool> &_undefs,
                   int _nCPUs, int _perm, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, _nCPUs, _perm, _last_seed), 
  CLUSTER_NOT_SIG(0),
  CLUSTER_HIGHHIGH(1),
  CLUSTER_LOWLOW(2),
  CLUSTER_OTHERPOS(3),
  CLUSTER_NEGATIVE(4),
  CLUSTER_UNDEFINED(5),
  CLUSTER_NEIGHBORLESS(6),
  data(_data), undefs(_undefs)
{
    labels.push_back("Not significant");
    labels.push_back("High-High");
    labels.push_back("Low-Low");
    labels.push_back("Other Positive");
    labels.push_back("Negative");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#b2182b");
    colors.push_back("#ef8a62");
    colors.push_back("#fddbc7");
    colors.push_back("#67adc7");
    colors.push_back("#464646");
    colors.push_back("#999999");

    GenUtils::StandardizeData(data, undefs);

    data_square.resize(num_obs, 0);
    for (int i=0; i<num_obs; ++i) {
        data_square[i] = data[i] * data[i];
    }

    Run();
}

UniGeary::~UniGeary() {

}

void UniGeary::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i]) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;

        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                // computer spatial lag
                double sp_lag = 0, sp_lag_square = 0;
                const std::vector<long>& nbrs = weights->GetNeighbors(i);
                unsigned int nn = 0;
                for (size_t j=0; j<nbrs.size(); ++j) {
                    if (nbrs[j] != i &&  !undefs[nbrs[j]]) { // not including the value at the location
                        sp_lag += data[ nbrs[j] ];
                        sp_lag_square += data_square[ nbrs[j] ];
                        nn += 1;
                    }
                }
                sp_lag = sp_lag / nn;
                sp_lag_square = sp_lag_square / nn;
                // compute geary's i
                lag_vec[i] = sp_lag;
                lisa_vec[i] = data_square[i] - 2.0 * data[i] * sp_lag + sp_lag_square;
                
                // assign the cluster
                if (data[i] > 0 && sp_lag > 0) cluster_vec[i] = CLUSTER_HIGHHIGH;
                else if (data[i] < 0 && sp_lag > 0) cluster_vec[i] = CLUSTER_OTHERPOS;
                else if (data[i] < 0 && sp_lag < 0) cluster_vec[i] = CLUSTER_LOWLOW;
                else cluster_vec[i] = CLUSTER_NEGATIVE; //data1[i] > 0 && sp_lag < 0
            }
        }
    }
}

void UniGeary::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA)
{
    int validNeighbors = 0;
    double permutedLag = 0;
    double permutedLagSquare = 0;
    int numNeighbors = permNeighbors.size();
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb]) {
            permutedLag += data[nb];
            permutedLagSquare += data_square[nb];
            validNeighbors ++;
        }
    }
    //NOTE: we shouldn't have to row-standardize or
    // multiply by data1[cnt]
    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors;
        permutedLagSquare /= validNeighbors;

    }

    const double localGearyPermuted = data_square[cnt] - 2.0 * data[cnt] * permutedLag + permutedLagSquare;
    permutedSA[perm] = localGearyPermuted;
}

uint64_t UniGeary::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
{
    double permGearySum = 0, permGearyMean = 0;
    for (int i=0; i<permutations; ++i) {
        permGearySum += permutedSA[i];
    }
    permGearyMean = permGearySum / permutations;
    uint64_t countLarger = 0;

    if (lisa_vec[cnt] <= permGearyMean) {
        // positive
        for (int i=0; i<permutations; ++i) {
            if (permutedSA[i] <= lisa_vec[cnt]) {
                countLarger += 1;
            }
            // positive HH cluster_vec[cnt] == 1CLUSTER_HIGHHIGH
            // positive LL cluster_vec[cnt] == 2CLUSTER_LOWLOW
            // positive && but in outlier qudrant: other pos
            if ((const unsigned long)cluster_vec[cnt] > CLUSTER_LOWLOW && (const unsigned long)cluster_vec[cnt] <
            CLUSTER_UNDEFINED) {
                cluster_vec[cnt] = CLUSTER_OTHERPOS;
            }
        }
    } else {
        // negative
        for (int i=0; i<permutations; ++i) {
            if (permutedSA[i] > lisa_vec[cnt]) {
                countLarger += 1;
            }
        }
        if ((const unsigned long)cluster_vec[cnt] < CLUSTER_UNDEFINED) {
            cluster_vec[cnt] = CLUSTER_NEGATIVE;
        }

    }
    return countLarger;
}

std::vector<int> UniGeary::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cuttoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if (sig_local_vec[i] > cuttoff &&
                (const unsigned long)cluster_vec[i] != CLUSTER_UNDEFINED &&
                (const unsigned long)cluster_vec[i] != CLUSTER_NEIGHBORLESS)
        {
            clusters[i] = CLUSTER_NOT_SIG;
        } else {
            clusters[i] = cluster_vec[i];
        }
    }
    return clusters;
}

