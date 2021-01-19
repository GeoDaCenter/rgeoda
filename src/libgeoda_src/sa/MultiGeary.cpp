//
// Created by Xun Li on 2019-12-01.
//

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"
#include "MultiGeary.h"

MultiGeary::MultiGeary(int num_obs, GeoDaWeight *w,
        const std::vector<std::vector<double> > &_data,
        const std::vector<std::vector<bool> > &_undefs,
        int _nCPUs, int _perm, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, _nCPUs, _perm, _last_seed),          
 CLUSTER_NOT_SIG(0),
 CLUSTER_POSITIVE(1),
 CLUSTER_NEGATIVE(2),
 CLUSTER_UNDEFINED(3),
 CLUSTER_NEIGHBORLESS(4),
 data(_data)
{
    labels.push_back("Not significant");
    labels.push_back("Positive");
    labels.push_back("Negative");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#336ea1");
    colors.push_back("#67adc7");
    colors.push_back("#464646");
    colors.push_back("#999999");

    std::vector<bool> undef_merge(num_obs, false);
    if (_undefs.size() > 0) {
        for (int i=0; i<num_obs; ++i) {
            for (size_t j = 0; j < _undefs.size(); ++j) {
                if ((int)_undefs[j].size() >= num_obs) {
                    break;
                }
                undef_merge[i] = undef_merge[i] || _undefs[j][i];
            }
        }
    }
    undefs = undef_merge;

    num_vars = data.size();
    for (int i=0; i < num_vars; ++i) {
        GenUtils::StandardizeData(data[i], undefs);
    }

    data_square.resize(num_vars);
    for (int i=0; i < num_vars; ++i) {
        data_square[i].resize(num_obs);
        for (int j = 0; j < num_obs; ++j) {
            data_square[i][j] = data[i][j] * data[i][j];
        }
    }

    Run();
}

MultiGeary::~MultiGeary() {

}

void MultiGeary::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i]) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;

        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                for (int v=0; v<num_vars; v++) {
                    // computer spatial lag
                    double sp_lag = 0, sp_lag_square = 0;
                    const std::vector<long>& nbrs = weights->GetNeighbors(i);
                    unsigned int nn = 0;

                    for (size_t j=0; j<nbrs.size(); ++j) {
                        if (nbrs[j] != i &&  !undefs[nbrs[j]]) { // not including the value at the location
                            sp_lag += data[v][ nbrs[j] ];
                            sp_lag_square += data_square[v][ nbrs[j] ];
                            nn += 1;
                        }
                    }
                    sp_lag = sp_lag / nn;
                    sp_lag_square = sp_lag_square / nn;
                    // compute geary's i
                    lag_vec[i] = sp_lag;
                    lisa_vec[i] += data_square[v][i] - 2.0 * data[v][i] * sp_lag + sp_lag_square;
                }

                lag_vec[i] /= num_vars;
                lisa_vec[i] /= num_vars;

                // assign the cluster not here but after permutation
            }
        }
    }
}

void MultiGeary::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA)
{
    int validNeighbors = 0;
    int numNeighbors = permNeighbors.size();

    // use permutation to compute the lag
    // compute the lag for binary weights
    std::vector<double> permutedLag(num_vars, 0);
    std::vector<double> permutedLagSquare(num_vars, 0);

    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb]) {
            validNeighbors ++;
            for (int v=0; v<num_vars; v++) {
                permutedLag[v] += data[v][nb];
                permutedLagSquare[v] += data_square[v][nb];
            }
        }
    }

    //NOTE: we shouldn't have to row-standardize or
    // multiply by data1[cnt]
    if (validNeighbors > 0 && row_standardize) {
        for (int v=0; v<num_vars; v++) {
            permutedLag[v] /= validNeighbors;
            permutedLagSquare[v] /= validNeighbors;
        }
    }

    double localGearyPermuted = 0;
    for (int v=0; v<num_vars; v++) {
        localGearyPermuted += data_square[v][cnt] - 2.0 * data[v][cnt] * permutedLag[v] + permutedLagSquare[v];
    }
    localGearyPermuted /= num_vars;
    permutedSA[perm] = localGearyPermuted;
}

uint64_t MultiGeary::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
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
            // ignore neighborless & undefined
            if ((const unsigned long)cluster_vec[cnt] < CLUSTER_UNDEFINED) {
                cluster_vec[cnt] = CLUSTER_POSITIVE;
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

std::vector<int> MultiGeary::GetClusterIndicators() {
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
