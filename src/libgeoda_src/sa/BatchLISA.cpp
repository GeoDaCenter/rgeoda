//
// Created by Xun Li on 2019-06-05.
//

#include <math.h>
#include <time.h>

#include "../GeoDaSet.h"
#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"

#include "BatchLISA.h"

#ifndef __USE_PTHREAD__
#include <boost/system/config.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#else
#include <pthread.h>

struct batchlisa_thread_args {
    BatchLISA *lisa;
    int start;
    int end;
    uint64_t seed_start;
};

void* batchlisa_thread_helper(void* voidArgs)
{
    batchlisa_thread_args *args = (batchlisa_thread_args*)voidArgs;
    args->lisa->CalcPseudoP_range(args->start, args->end, args->seed_start);
    return 0;
}
#endif


BatchLISA::BatchLISA(int num_obs, GeoDaWeight* w,
                     const std::vector<std::vector<bool> >& _undefs,
                     double _significance_cutoff,
                     int _nCPUs, int _perm, uint64_t _last_seed)
    :  nCPUs(_nCPUs),
    num_obs(num_obs),
    row_standardize(true),
    permutations(_perm),
    significance_cutoff(_significance_cutoff),
    user_sig_cutoff(0),
    has_undefined(false),
    has_isolates(w->HasIsolates()),
    calc_significances(true),
    last_seed_used(_last_seed),
    reuse_last_seed(true),
    weights(w),
    undefs(_undefs)
{
    //SetSignificanceFilter(1);
}

BatchLISA::~BatchLISA()
{
}

void BatchLISA::Run()
{
    sig_local_vec.resize(num_batch);
    sig_cat_vec.resize(num_batch);
    cluster_vec.resize(num_batch);
    lag_vec.resize(num_batch);
    lisa_vec.resize(num_batch);

    for (int i=0; i<num_batch; ++i) {
        sig_local_vec[i].resize(num_obs, 0);
        sig_cat_vec[i].resize(num_obs, 0);
        cluster_vec[i].resize(num_obs, 0);
        lag_vec[i].resize(num_obs, 0);
        lisa_vec[i].resize(num_obs, 0);
    }

    nn_vec.resize(num_obs, 0);
    for (int i=0; i<num_obs; i++) {
        nn_vec[i] = weights->GetNbrSize(i);
    }

    ComputeLoalSA();

    if (calc_significances) {
        CalcPseudoP();
    }
}

void BatchLISA::SetSignificanceFilter(int filter_id)
{
    if (filter_id == -1) {
        // user input cutoff
        significance_filter = filter_id;
        return;
    }
    // 0: >0.05 1: 0.05, 2: 0.01, 3: 0.001, 4: 0.0001
    if (filter_id < 1 || filter_id > 4) return;
    significance_filter = filter_id;
    if (filter_id == 1) significance_cutoff = 0.05;
    if (filter_id == 2) significance_cutoff = 0.01;
    if (filter_id == 3) significance_cutoff = 0.001;
    if (filter_id == 4) significance_cutoff = 0.0001;
}

int BatchLISA::GetSignificanceFilter()
{
    return significance_filter;
}

double BatchLISA::GetSignificanceCutoff()
{
    return significance_cutoff;
}

void BatchLISA::SetSignificanceCutoff(double val)
{
    significance_cutoff = val;
}

double BatchLISA::GetUserCutoff()
{
    return user_sig_cutoff;
}

void BatchLISA::SetUserCutoff(double val)
{
    user_sig_cutoff = val;
}

double BatchLISA::GetFDR(double current_p, int idx)
{
    // bound check
    if (idx < 0 || idx >= (int)sig_local_vec.size()-1)
        return 0;

    double fdr = 0; //False Discovery Rate
    std::vector<double> pvals = sig_local_vec[idx]; // make sure copy
    // FDR
    // sort all p-values from smallest to largets
    std::sort(pvals.begin(), pvals.end());

    int i_0 = -1;
    bool stop = false;
    double p_start = current_p;

    while (!stop) {
        // find the i_0 that corresponds to p = alpha
        for (int i=1; i<=num_obs; i++) {
            if (pvals[i] >= p_start) {
                if (i_0 == i) {
                    stop = true;
                }
                i_0 = i;
                break;
            }
        }
        if (i_0 < 0)
            stop = true;

        // compute p* = i_0 x alpha / N
        p_start = i_0 * current_p / (double)num_obs ;
    }

    if (i_0 < 0)
        p_start = 0.0;

    fdr = p_start;

    return fdr;
}


double BatchLISA::GetBO(double current_p)
{
    //double bo; //Bonferroni bound
    double bonferroni_bound = current_p / (double)num_obs;
    return bonferroni_bound;
}

int BatchLISA::GetNumPermutations()
{
    return permutations;
}

void BatchLISA::SetNumPermutations(int val)
{
    permutations = val;
}

uint64_t BatchLISA::GetLastUsedSeed()
{
    return last_seed_used;
}

bool BatchLISA::IsReuseLastSeed()
{
    return reuse_last_seed;
}

void BatchLISA::SetReuseLastSeed(bool reuse)
{
    reuse_last_seed = reuse;
}

void BatchLISA::SetLastUsedSeed(uint64_t seed)
{
    reuse_last_seed = true;
    last_seed_used = seed;
}

bool BatchLISA::GetHasIsolates()
{
    return has_isolates;
}

bool BatchLISA::GetHasUndefined()
{
    return has_undefined;

}

void BatchLISA::CalcPseudoP()
{
    if (!calc_significances) return;
    CalcPseudoP_threaded();
}

void BatchLISA::CalcPseudoP_threaded()
{
#ifndef __USE_PTHREAD__
    if (nCPUs <= 0) nCPUs = boost::thread::hardware_concurrency();
    boost::thread_group threadPool;
#else
    pthread_t* threadPool = new pthread_t[nCPUs];
    struct batchlisa_thread_args* args = new batchlisa_thread_args[nCPUs];
#endif

    // divide up work according to number of observations
    // and number of CPUs
    int work_chunk = num_obs / nCPUs;

    if (work_chunk == 0) work_chunk = 1;

    //int obs_start = 0;
    //int obs_end = obs_start + work_chunk;

    int quotient = num_obs / nCPUs;
    int remainder = num_obs % nCPUs;
    int tot_threads = (quotient > 0) ? nCPUs : remainder;

    if (!reuse_last_seed) last_seed_used = time(0);

    for (int i=0; i<tot_threads; i++) {
        int a=0;
        int b=0;
        if (i < remainder) {
            a = i*(quotient+1);
            b = a+quotient;
        } else {
            a = remainder*(quotient+1) + (i-remainder)*quotient;
            b = a+quotient-1;
        }
        uint64_t seed_start = last_seed_used+a;
        //uint64_t seed_end = seed_start + ((uint64_t) (b-a));

#ifndef __USE_PTHREAD__
        boost::thread* worker = new boost::thread(boost::bind(&BatchLISA::CalcPseudoP_range,this, a, b, seed_start));
        threadPool.add_thread(worker);
#else
        args[i].lisa = this;
        args[i].start = a;
        args[i].end = b;
        args[i].seed_start = seed_start;
        if (pthread_create(&threadPool[i], NULL, &batchlisa_thread_helper, &args[i])) {
            perror("Thread create failed.");
        }
#endif
    }
#ifndef __USE_PTHREAD__
    threadPool.join_all();
#else
    for (int j = 0; j < nCPUs; j++) {
        pthread_join(threadPool[j], NULL);
    }
    delete[] args;
    delete[] threadPool;
#endif
}

void BatchLISA::CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    GeoDaSet workPermutation(num_obs);
    int max_rand = num_obs-1;
    //bool perm_valid = true;
    int numNeighbors;

    for (int cnt=obs_start; cnt<=obs_end; cnt++) {
        numNeighbors = weights->GetNbrSize(cnt);
        if (numNeighbors == 0) {
            for (int  v=0; v < num_batch; ++v) {
                sig_cat_vec[v][cnt] = 5; // neighborless cat
            }
        } else {
            std::vector<std::vector<double> > permutedSA(num_batch);
            for (int i=0; i<num_batch; ++i) permutedSA[i].resize(permutations);

            for (int perm = 0; perm < permutations; perm++) {
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

                PermLocalSA(cnt, perm, permNeighbors, permutedSA);
            }

            std::vector<uint64_t> countLarger = CountLargerSA(cnt, permutedSA);

            for (int v=0; v < num_batch; ++v) {
                double _sigLocal = (countLarger[v] + 1.0) / (permutations + 1);

                // 'significance' of local sa
                if (_sigLocal <= 0.0001) sig_cat_vec[v][cnt] = 4;
                else if (_sigLocal <= 0.001) sig_cat_vec[v][cnt] = 3;
                else if (_sigLocal <= 0.01) sig_cat_vec[v][cnt] = 2;
                else if (_sigLocal <= 0.05) sig_cat_vec[v][cnt] = 1;
                else sig_cat_vec[v][cnt] = 0;

                if (undefs[v][cnt]) sig_cat_vec[v][cnt] = 6; // undefined

                sig_local_vec[v][cnt] = _sigLocal;
                // observations with no neighbors get marked as isolates
                // NOTE: undefined should be marked as well, however, since undefined_cat has covered undefined category,
                // we don't need to handle here
            }
        }
    }
}


std::vector<std::string> BatchLISA::GetDefaultCategories()
{
    std::vector<std::string> cats;
    cats.push_back("p = 0.05");
    cats.push_back("p = 0.01");
    cats.push_back("p = 0.001");
    cats.push_back("p = 0.0001");
    return cats;
}

std::vector<double> BatchLISA::GetDefaultCutoffs()
{
    std::vector<double> cutoffs;
    cutoffs.push_back(0.05);
    cutoffs.push_back(0.01);
    cutoffs.push_back(0.001);
    cutoffs.push_back(0.0001);
    return cutoffs;
}

std::vector<double> BatchLISA::GetLocalSignificanceValues(int idx)
{
    return sig_local_vec[idx];
}

std::vector<int> BatchLISA::GetClusterIndicators(int idx)
{
    return cluster_vec[idx];
}

std::vector<int> BatchLISA::GetSigCatIndicators(int idx)
{
    return sig_cat_vec[idx];
}

std::vector<int> BatchLISA::GetNumNeighbors()
{
    return nn_vec;
}

std::vector<double> BatchLISA::GetSpatialLagValues(int idx)
{
    return lag_vec[idx];
}

std::vector<double> BatchLISA::GetLISAValues(int idx)
{
    return lisa_vec[idx];
}

bool BatchLISA::IsRowStandardize() const {
    return row_standardize;
}

void BatchLISA::SetRowStandardize(bool rowStandardize) {
    row_standardize = rowStandardize;
}

int BatchLISA::GetNumThreads() const {
    return nCPUs;
}

void BatchLISA::SetNumThreads(int n_threads) {
    nCPUs = n_threads;
}

std::vector<std::string> BatchLISA::GetLabels()
{
    return labels;
}

std::vector<std::string> BatchLISA::GetColors()
{
    return colors;
}
