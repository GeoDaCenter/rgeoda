//
// Created by Xun Li on 2019-12-01.
//

#include <math.h>
#include <vector>

#include "../GeoDaSet.h"
#include "../weights/GeodaWeight.h"
#include "../GenUtils.h"
#include "MultiJoinCount.h"

MultiJoinCount::MultiJoinCount(int num_obs, GeoDaWeight *w,
                           const std::vector<std::vector<double> > &_data,
                           const std::vector<std::vector<bool> > &_undefs,
                               int _nCPUs, int _perm, uint64_t _last_seed)
        : LISA(num_obs, w, _undefs, _nCPUs, _perm, _last_seed), 
          CLUSTER_NOT_SIG(0),
          CLUSTER_SIG(1),
          CLUSTER_UNDEFINED(2),
          CLUSTER_NEIGHBORLESS(3),
          data(_data)
{
    labels.push_back("Not significant");
    labels.push_back("Significant");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#348124"); // green
    colors.push_back("#464646");
    colors.push_back("#999999");

    num_vars = data.size();

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

    zz.resize(num_obs, 1);
    for (int i=0; i<num_obs; i++) {
        for (int v = 0; v < num_vars; v++) {
            zz[i] = zz[i] * data[v][i];  // 0 or 1
        }
    }

    Run();
}

MultiJoinCount::~MultiJoinCount() {

}

void MultiJoinCount::ComputeLoalSA() {

    int sum = 0;
    for (int i=0; i<num_obs; i++) {
        if (!undefs[i])
            sum += zz[i];
    }
    bool nocolocation = sum == 0;

    // bivariate local join count -- colocation and no-colocation
    // multivariate local join count -- colocation only
    if (nocolocation) {
        // here only bivariate apply to no-colocation case
        // why? think about if add a third variable, it would break
        // the no-colocation case!
        for (int i=0; i<num_obs; i++) {
            if (undefs[i] == true) {
                zz[i] = 0;
                continue;
            }
            zz[i] = data[1][i];
        }
        for (int i=0; i<num_obs; i++) {
            if (undefs[i]) {
                lag_vec[i] = 0;
                lisa_vec[i] = 0;
                cluster_vec[i] = CLUSTER_UNDEFINED;
            } else {
                if (weights->GetNbrSize(i) == 0) {
                    cluster_vec[i] = CLUSTER_NEIGHBORLESS;
                } else {
                    if (data[0][i] > 0) { // x_i = 1
                        int nbr_size = weights->GetNbrSize(i);
                        const std::vector<long>& nbrs = weights->GetNeighbors(i);
                        for (int j=0; j<nbr_size; ++j) {
                            if (nbrs[j] != i &&  !undefs[nbrs[j]])
                                lisa_vec[i] += zz[ nbrs[j] ];
                        }
                    }
                }
            }
        }
    } else {
        // could be both bivariate and multivariate cases
        for (int i=0; i<num_obs; i++) {
            if (undefs[i]) {
                lag_vec[i] = 0;
                lisa_vec[i] = 0;
                cluster_vec[i] = CLUSTER_UNDEFINED;
            } else {
                if (zz[i] > 0) { // x_i.z_i = 1
                    int nbr_size = weights->GetNbrSize(i);
                    const std::vector<long> &nbrs = weights->GetNeighbors(i);
                    for (int j = 0; j < nbr_size; ++j) {
                        if (nbrs[j] != i && !undefs[nbrs[j]]) {
                            // compute the number of neighbors with
                            // x_j.z_j = 1 (zz=1) as a spatial lag
                            lisa_vec[i] += zz[nbrs[j]];
                        }
                    }
                }
            }
        }
    }

}

void MultiJoinCount::CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    GeoDaSet workPermutation(num_obs);
    int max_rand = num_obs-1;

    for (int cnt=obs_start; cnt<=obs_end; cnt++) {
        if (undefs[cnt]) {
            sig_cat_vec[cnt] = 6; // undefined cat
            continue;
        }
        if (lisa_vec[cnt] == 0) {
            sig_local_vec[cnt] = 0;
            continue;
        }

        // get full neighbors even if has undefined value
        int numNeighbors = weights->GetNbrSize(cnt);
        if (numNeighbors == 0) {
            sig_cat_vec[cnt] = 5; // neighborless cat
            // isolate: don't do permutation
            continue;
        }

        std::vector<double> permutedSA(permutations, 0);
        for (int perm=0; perm<permutations; perm++) {
            int rand=0, newRandom;
            double rng_val;
            while (rand < numNeighbors) {
                // computing 'perfect' permutation of given size
                rng_val = Gda::ThomasWangHashDouble(seed_start++) * max_rand;
                // round is needed to fix issue
                // https://github.com/GeoDaCenter/geoda/issues/488
                newRandom = (int)(rng_val<0.0?ceil(rng_val - 0.5):floor(rng_val + 0.5));

                if (newRandom != cnt && !workPermutation.Belongs(newRandom) && weights->GetNbrSize(newRandom)>0) {
                    workPermutation.Push(newRandom);
                    rand++;
                }
            }
            std::vector<int> permNeighbors(numNeighbors);
            for (int cp=0; cp<numNeighbors; cp++) {
                permNeighbors[cp] = workPermutation.Pop();
            }

            PermLocalSA(cnt, perm, permNeighbors, permutedSA);
        }

        uint64_t countLarger = CountLargerSA(cnt, permutedSA);
        double _sigLocal = (countLarger+1.0)/(permutations+1);

        // 'significance' of local sa
        if (_sigLocal <= 0.0001) sig_cat_vec[cnt] = 4;
        else if (_sigLocal <= 0.001) sig_cat_vec[cnt] = 3;
        else if (_sigLocal <= 0.01) sig_cat_vec[cnt] = 2;
        else if (_sigLocal <= 0.05) sig_cat_vec[cnt] = 1;
        else sig_cat_vec[cnt] = 0;

        sig_local_vec[cnt] = _sigLocal;
        // observations with no neighbors get marked as isolates
        // NOTE: undefined should be marked as well, however, since undefined_cat has covered undefined category, we don't need to handle here
    }
}

void
MultiJoinCount::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double> &permutedSA) {

    int validNeighbors = 0;
    double permutedLag = 0;
    int numNeighbors = permNeighbors.size();
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb]) {
            permutedLag += zz[nb];
            validNeighbors ++;
        }
    }
    permutedSA[perm] = permutedLag;
}

uint64_t MultiJoinCount::CountLargerSA(int cnt, const std::vector<double> &permutedSA) {
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

std::vector<int> MultiJoinCount::GetClusterIndicators() {
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
