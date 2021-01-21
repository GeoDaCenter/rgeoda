#include <iostream>
#include <float.h>

#include "weights/GeodaWeight.h"
#include "sa/UniGeary.h"
#include "sa/UniG.h"
#include "sa/UniGstar.h"
#include "sa/UniJoinCount.h"
#include "sa/UniLocalMoran.h"
#include "sa/MultiGeary.h"
#include "sa/MultiJoinCount.h"
#include "sa/BatchLocalMoran.h"
#include "GenUtils.h"
#include "gda_sa.h"

LISA *gda_localg(GeoDaWeight *w,
                 const std::vector<double> &data,
                 const std::vector<bool> &undefs,
                 double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniG *localg = new UniG(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    localg->SetSignificanceCutoff(significance_cutoff);
    return localg;
}

LISA *gda_localgstar(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniGstar *localgstar = new UniGstar(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    localgstar->SetSignificanceCutoff(significance_cutoff);
    return localgstar;
}

LISA *gda_localmoran(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniLocalMoran *lisa = new UniLocalMoran(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    lisa->SetSignificanceCutoff(significance_cutoff);
    return lisa;
}

BatchLISA *gda_batchlocalmoran(GeoDaWeight *w,
                               const std::vector<std::vector<double> > &data,
                               const std::vector<std::vector<bool> > &undefs,
                               double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<std::vector<bool> > copy_undefs = undefs;
    if (undefs.empty())
    {
        copy_undefs.resize(data.size());
        for (size_t i = 0; i < data.size(); ++i)
            copy_undefs[i].resize(num_obs, false);
    }

    BatchLISA *bm = new BatchLocalMoran(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    bm->SetSignificanceCutoff(significance_cutoff);
    return bm;
}

LISA *gda_localgeary(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniGeary *geary = new UniGeary(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    geary->SetSignificanceCutoff(significance_cutoff);
    return geary;
}

LISA *gda_localmultigeary(GeoDaWeight *w,
                     const std::vector<std::vector<double> > &data,
                     const std::vector<std::vector<bool> > &undefs,
                     double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    MultiGeary *geary = new MultiGeary(num_obs, w, data, undefs, nCPUs, perm, last_seed);
    geary->SetSignificanceCutoff(significance_cutoff);
    return geary;
}

LISA *gda_localjoincount(GeoDaWeight *w,
                    const std::vector<double> &data,
                    const std::vector<bool> &undefs,
                    double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniJoinCount *jc = new UniJoinCount(num_obs, w, data, copy_undefs, nCPUs, perm, last_seed);
    jc->SetSignificanceCutoff(significance_cutoff);
    return jc;
}

LISA *gda_localmultijoincount(GeoDaWeight *w,
                         const std::vector<std::vector<double> > &data,
                         const std::vector<std::vector<bool> > &undefs,
                         double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    MultiJoinCount *jc = new MultiJoinCount(num_obs, w, data, undefs, nCPUs, perm, last_seed);
    jc->SetSignificanceCutoff(significance_cutoff);
    return jc;
}

double gda_fdr(LISA *lisa, double current_p)
{
    if (lisa == 0)
        return 0;

    return lisa->GetFDR(current_p);
}

double gda_bo(LISA *lisa, double current_p)
{
    if (lisa == 0)
        return 0;

    return lisa->GetBO(current_p);
}

LISA *gda_quantilelisa(GeoDaWeight *w, unsigned int k, unsigned int quantile, const std::vector<double> &data,
                       const std::vector<bool> &undefs,
                       double significance_cutoff, int nCPUs, int perm, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    if (k < 1 || (int)k >= num_obs)
        return 0;

    if (quantile < 0 || quantile > k)
        return 0;

    std::vector<bool> copy_undefs = undefs; //copy
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }

    std::vector<double> breaks = GenUtils::QuantileBreaks(k, data, copy_undefs);

    quantile = quantile - 1;
    double break_left = DBL_MIN;
    double break_right = DBL_MAX;

    if (quantile == 0)
    {
        break_right = breaks[quantile];
    }
    else if (quantile == breaks.size())
    {
        break_left = breaks[quantile - 1];
    }
    else
    {
        break_left = breaks[quantile - 1];
        break_right = breaks[quantile];
    }

    std::vector<double> bin_data(num_obs, 0);

    for (int i = 0; i < num_obs; ++i)
    {
        if (data[i] >= break_left && data[i] < break_right)
        {
            bin_data[i] = 1;
        }
    }

    // apply local join count on binary data
    UniJoinCount *jc = new UniJoinCount(num_obs, w, bin_data, copy_undefs, nCPUs, perm, last_seed);
    jc->SetSignificanceCutoff(significance_cutoff);
    return jc;
}

LISA *gda_multiquantilelisa(GeoDaWeight *w, const std::vector<int>& k_s, const std::vector<int>& quantile_s, const std::vector<std::vector<double> > &data_s,
                       const std::vector<std::vector<bool> > &undefs_s, double significance_cutoff,
                       int nCPUs, int permutations, int last_seed_used)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    if (k_s.size() != quantile_s.size() || k_s.size() != data_s.size()) 
        return 0;

    // multi local joincount
    int num_vars = k_s.size();
    std::vector<std::vector<double> > data; // for storing binary data
    std::vector<std::vector<bool> > undefs = undefs_s;
    if (undefs.empty()) {
        undefs.resize(num_vars);
        for (size_t i=0; i< undefs.size(); ++i) {
            undefs[i].resize(num_obs, false);
        }
    }

    for (int i=0; i<num_vars; ++i) {
        int k = k_s[i];
        int q = quantile_s[i];

        std::vector<bool> copy_undefs = undefs[i];
        std::vector<double> copy_data = data_s[i];

        std::vector<double> breaks = GenUtils::QuantileBreaks(k, copy_data, copy_undefs); 

        q = q - 1;
        double break_left = DBL_MIN;
        double break_right = DBL_MAX;

        if (q == 0)
        {
            break_right = breaks[q];
        }
        else if (q== breaks.size())
        {
            break_left = breaks[q - 1];
        }
        else
        {
            break_left = breaks[q - 1];
            break_right = breaks[q];
        }

        std::vector<double> bin_data(num_obs, 0);

        for (int j = 0; j < num_obs; ++j)
        {
            if (copy_data[j] >= break_left && copy_data[j] < break_right)
            {
                bin_data[j] = 1;
            }
        }

        data.push_back(bin_data);
    }

    MultiJoinCount *jc = new MultiJoinCount(num_obs, w, data, undefs, nCPUs, permutations, last_seed_used);
    jc->SetSignificanceCutoff(significance_cutoff);

    return jc;
}
