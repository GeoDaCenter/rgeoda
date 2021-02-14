//
// Created by Xun Li on 9/27/19.
//
#include <vector>
#include <iostream>
#include <math.h>

#include "UniJoinCount.h"
#include "../weights/GeodaWeight.h"
#include "../weights/GwtWeight.h"
#include "../GenUtils.h"
#include "../GeoDaSet.h"

UniJoinCount::UniJoinCount(int num_obs, GeoDaWeight *w,
        const std::vector<double> &_data,
        const std::vector<bool> &_undefs,
        double significance_cutoff,
        int _nCPUs, int _perm, const std::string& _permutation_method, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, significance_cutoff, _nCPUs, _perm, _permutation_method, _last_seed),
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
                undefs[i] = true; // same logic as in GeoDa
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

void UniJoinCount::CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    GeoDaSet workPermutation(num_obs);
    int max_rand = num_obs-1;

#ifdef __JSGEODA__
    std::string wuid = weights->GetUID();
    bool using_cache = has_cached_perm[wuid];
    std::vector<std::vector<int> >& cache = cached_perm_nbrs[wuid];
#endif

    for (int cnt=obs_start; cnt<=obs_end; cnt++) {
        if (undefs[cnt]) {
            sig_cat_vec[cnt] = 6; // undefined
            continue;
        }

        if (lisa_vec[cnt] == 0) {
            sig_local_vec[cnt] = -1.0;
            continue;
        }

        // get full neighbors even if has undefined value
        int numNeighbors = weights->GetNbrSize(cnt);
        if (numNeighbors <= 0) {
            sig_cat_vec[cnt] = 5; // neighborless cat
            // isolate: don't do permutation
            continue;
        }

#ifdef __JSGEODA__
        if (using_cache == false) {
            for (size_t perm = 0; perm < permutations; perm++) {
                int rand = 0, newRandom;
                double rng_val;
                while (rand < numNeighbors) {
                    // computing 'perfect' permutation of given size
                    rng_val = Gda::ThomasWangHashDouble(seed_start++) * max_rand;
                    // round is needed to fix issue
                    // https://github.com/GeoDaCenter/geoda/issues/488
                    newRandom = (int) (rng_val < 0.0 ? ceil(rng_val - 0.5) : floor(rng_val + 0.5));

                    if (newRandom != cnt && !workPermutation.Belongs(newRandom) && weights->GetNbrSize(newRandom) > 0) {
                        workPermutation.Push(newRandom);
                        rand++;
                    }
                }
                std::vector<int> permNeighbors(numNeighbors);
                for (int cp = 0; cp < numNeighbors; cp++) {
                    permNeighbors[cp] = workPermutation.Pop();
                }
                cache.push_back(permNeighbors);
                PermLocalSA(cnt, perm, permNeighbors, permutedSA);

            }
        } else {
            for (size_t perm = 0; perm < permutations; perm++) {
                PermLocalSA(cnt, perm, cache[perm], permutedSA);
            }
        }
#else

        int countLarger = 0;
        for (int perm=0; perm<permutations; perm++) {
            int rand=0, newRandom;
            double rng_val;
            while (rand < numNeighbors) {
                // computing 'perfect' permutation of given size
                rng_val = Gda::ThomasWangHashDouble(seed_start++) * max_rand;
                // round is needed to fix issue
                // https://github.com/GeoDaCenter/geoda/issues/488
                newRandom = (int)(rng_val<0.0?ceil(rng_val - 0.5):floor(rng_val + 0.5));

                if (newRandom != cnt && !workPermutation.Belongs(newRandom) && undefs[newRandom] == false) {
                    workPermutation.Push(newRandom);
                    rand++;
                }
            }

            double perm_jc = 0;
            // use permutation to compute the lags
            for (int cp=0; cp<numNeighbors; cp++) {
                int perm_idx = workPermutation.Pop();
                perm_jc += data[perm_idx];
            }

            // binary weights
            if (perm_jc >= lisa_vec[cnt]) {
                countLarger++;
            }
        }
#endif

        // pick the smallest counts
        if (permutations-countLarger <= countLarger) {
            countLarger = permutations-countLarger;
        }

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


void UniJoinCount::PermCalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    for (int cnt=obs_start; cnt<=obs_end; cnt++) {
        if (undefs[cnt]) {
            sig_cat_vec[cnt] = 6; // undefined
            continue;
        }
        // ignore local join count == 0
        if (lisa_vec[cnt] == 0) {
            sig_local_vec[cnt] = -1.0;
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
        for (size_t perm = 0; perm < permutations; perm++) {
            PermLocalSA(cnt, perm, numNeighbors, perm_table[perm], permutedSA);
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
    }
}

void UniJoinCount::PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors,
                                std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (nb >= cnt) nb = nb + 1; // self "cnt" should be excluded, index should be adjusted
        if (!undefs[nb]) {
            permutedLag += data[nb];
            validNeighbors ++;
        }
    }
    permutedSA[perm] = permutedLag;
}

void UniJoinCount::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>
        &permutedSA) {
    // not used, see  UniJoinCount::CalcPseudoP_range
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
    double cutoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if (sig_local_vec[i] <= cutoff ) {
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

