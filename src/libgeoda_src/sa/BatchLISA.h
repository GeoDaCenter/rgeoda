//
// Created by Xun Li on 2019-06-05.
//

#ifndef GEODA_BATCHLISA_H
#define GEODA_BATCHLISA_H

#ifdef _WIN32
// for uint64_t typedef
#include "../pg/geoms.h"
#endif

#include <list>
#include <string>

class GeoDaWeight;

class BatchLISA
{
public:
    BatchLISA(int num_obs, GeoDaWeight* w, const std::vector<std::vector<bool> >& undefs,
            double significance_cutoff,
         int nCPUs,
         int permutations,
         uint64_t last_seed_used);

    virtual ~BatchLISA();

    virtual void ComputeLoalSA() = 0;

    virtual void CalcPseudoP();

    virtual void CalcPseudoP_threaded();

    virtual void CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start);

    // compare local SA value of current obs to local SA values of random picked nbrs
    virtual void PermLocalSA(int cnt, int perm,
            const std::vector<int>& permNeighbors,
            std::vector<std::vector<double> >& permutedSA) = 0;

    virtual std::vector<uint64_t> CountLargerSA(int cnt,
            const std::vector<std::vector<double> >& permutedSA) = 0;

    virtual void Run();

    virtual void SetSignificanceFilter(int filter_id);
    virtual int GetSignificanceFilter();

    virtual double GetSignificanceCutoff();
    virtual void SetSignificanceCutoff(double val);

    virtual double GetUserCutoff();
    virtual void SetUserCutoff(double val);

    virtual double GetBO(double current_p);

    virtual double GetFDR(double current_p, int idx);

    virtual int GetNumPermutations();
    virtual void SetNumPermutations(int val);

    virtual uint64_t GetLastUsedSeed();
    virtual void SetLastUsedSeed(uint64_t seed);

    virtual bool IsReuseLastSeed();
    virtual void SetReuseLastSeed(bool reuse);

    virtual bool GetHasIsolates();

    virtual bool GetHasUndefined();

    virtual std::vector<std::string> GetDefaultCategories();

    virtual std::vector<double> GetDefaultCutoffs();

    virtual std::vector<double> GetLocalSignificanceValues(int idx);

    virtual std::vector<int> GetClusterIndicators(int idx);

    virtual std::vector<int> GetSigCatIndicators(int idx);

    virtual std::vector<int> GetNumNeighbors();

    virtual std::vector<double> GetSpatialLagValues(int idx);

    virtual std::vector<double> GetLISAValues(int idx);

    virtual bool IsRowStandardize() const;

    virtual void SetRowStandardize(bool rowStandardize);

    virtual int GetNumThreads() const;

    virtual void SetNumThreads(int n_threads);

    virtual std::vector<std::string> GetLabels();

    virtual std::vector<std::string> GetColors();

protected:
    int nCPUs;

    int num_obs; // total # obs including neighborless obs

    bool row_standardize;

    int significance_filter; // 0: >0.05 1: 0.05, 2: 0.01, 3: 0.001, 4: 0.0001

    int permutations; // any number from 9 to 99999, 99 will be default

    double significance_cutoff; // either 0.05, 0.01, 0.001 or 0.0001

    double user_sig_cutoff; // user defined cutoff

    bool has_undefined;

    bool has_isolates;

    bool calc_significances; // if false, then p-vals will never be needed

    uint64_t last_seed_used;

    bool reuse_last_seed;

    GeoDaWeight* weights;

    int num_batch;

    // return values for all input variables
    std::vector<std::vector<bool> > undefs;
    std::vector<std::vector<double> > sig_local_vec;
    std::vector<std::vector<int> > sig_cat_vec;
    std::vector<std::vector<int> > cluster_vec;
    std::vector<std::vector<double> > lag_vec;
    std::vector<std::vector<double> > lisa_vec;

    std::vector<int> nn_vec;
    std::vector<std::string> labels;
    std::vector<std::string> colors;
};


#endif //GEODA_BATCHLISA_H
