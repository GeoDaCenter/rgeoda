#include <iostream>
#include <float.h>
#include <boost/algorithm/string.hpp>

#include "knn/ANN/ANN.h"
#include "weights/GeodaWeight.h"
#include "weights/GalWeight.h"
#include "sa/UniGeary.h"
#include "sa/UniG.h"
#include "sa/UniGstar.h"
#include "sa/UniJoinCount.h"
#include "sa/UniLocalMoran.h"
#include "sa/MultiGeary.h"
#include "sa/MultiJoinCount.h"
#include "sa/BatchLocalMoran.h"
#include "gda_data.h"
#include "GenUtils.h"
#include "gda_sa.h"
#include "gda_weights.h"
#include "gda_interface.h"

LISA *gda_localg(GeoDaWeight *w,
                 const std::vector<double> &data,
                 const std::vector<bool> &undefs,
                 double significance_cutoff, int nCPUs, int perm,
                 const std::string& perm_method, int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniG *localg = new UniG(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return localg;
}

LISA *gda_localgstar(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniGstar *localgstar = new UniGstar(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return localgstar;
}

LISA *gda_localmoran(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniLocalMoran *lisa = new UniLocalMoran(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return lisa;
}

LISA *gda_localmoran_eb(GeoDaWeight *w,
                        const std::vector<double> &event_data,
                        const std::vector<double> &base_data,
                        double significance_cutoff,
                        int nCPUs, int permutations, const std::string& perm_method, int last_seed_used)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<double> E, P, local_smoothed_results(num_obs, 0.0);
    E = event_data;
    P = base_data;

    // eb_rate_standardized
    std::vector<bool> undef_res(num_obs, false);
    bool success = gda_rateStandardizeEB(P, E, local_smoothed_results, undef_res);

    if (success == false)
        return 0;

    UniLocalMoran *lisa = new UniLocalMoran(num_obs, w, local_smoothed_results, undef_res, significance_cutoff, nCPUs, permutations, perm_method, last_seed_used);
    return lisa;
}

BatchLISA *gda_batchlocalmoran(GeoDaWeight *w,
                               const std::vector<std::vector<double> > &data,
                               const std::vector<std::vector<bool> > &undefs,
                               double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
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

    BatchLISA *bm = new BatchLocalMoran(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, last_seed);
    return bm;
}

LISA *gda_localgeary(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs,
                     double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniGeary *geary = new UniGeary(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return geary;
}

LISA *gda_localmultigeary(GeoDaWeight *w,
                     const std::vector<std::vector<double> > &data,
                     const std::vector<std::vector<bool> > &undefs,
                     double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    MultiGeary *geary = new MultiGeary(num_obs, w, data, undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return geary;
}

LISA *gda_localjoincount(GeoDaWeight *w,
                    const std::vector<double> &data,
                    const std::vector<bool> &undefs,
                    double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty())
    {
        copy_undefs.resize(num_obs, false);
    }
    UniJoinCount *jc = new UniJoinCount(num_obs, w, data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return jc;
}

LISA *gda_localmultijoincount(GeoDaWeight *w,
                         const std::vector<std::vector<double> > &data,
                         const std::vector<std::vector<bool> > &undefs,
                         double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
{
    if (w == 0)
        return 0;

    int num_obs = w->num_obs;

    MultiJoinCount *jc = new MultiJoinCount(num_obs, w, data, undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
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
                       double significance_cutoff, int nCPUs, int perm, const std::string& perm_method,  int last_seed)
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
    UniJoinCount *jc = new UniJoinCount(num_obs, w, bin_data, copy_undefs, significance_cutoff, nCPUs, perm, perm_method, last_seed);
    return jc;
}

LISA *gda_multiquantilelisa(GeoDaWeight *w, const std::vector<int>& k_s, const std::vector<int>& quantile_s, const std::vector<std::vector<double> > &data_s,
                       const std::vector<std::vector<bool> > &undefs_s, double significance_cutoff,
                       int nCPUs, int permutations, const std::string& perm_method, int last_seed_used)
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
        else if (q == (int)breaks.size())
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

    MultiJoinCount *jc = new MultiJoinCount(num_obs, w, data, undefs, significance_cutoff, nCPUs, permutations, perm_method, last_seed_used);
    return jc;
}

std::vector<std::vector<double> > gda_neighbor_match_test(AbstractGeoDa* geoda, unsigned int knn,
                                                          double power,
                                                          bool is_inverse,
                                                          bool is_arc,
                                                          bool is_mile,
                                                          const std::vector<std::vector<double> >& _data,
                                                          const std::string& scale_method,
                                                          const std::string& dist_type)
{
    int rows = geoda->GetNumObs();
    int columns = _data.size();

    // knn weights
    std::string kernel = "";
    double bandwidth = 0;
    bool adaptive_bandwidth = false;
    bool use_kernel_diagonal = false;
    std::string polyid = "";
    GeoDaWeight* sw = gda_knn_weights(geoda, knn, power, is_inverse, is_arc,
                                      is_mile, kernel, bandwidth, adaptive_bandwidth,
                                      use_kernel_diagonal, polyid);

    // transform data
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    double **input_data = new double*[rows];
    for (int i=0; i<rows; i++) {
        input_data[i] = new double[columns];
    }
    for (int i=0; i<columns; ++i) {
        for (int k=0; k< rows;k++) { // row
            input_data[k][i] = data[i][k];
        }
    }


    // create knn variable weights
    double eps = 0; // error bound
    ANN_DIST_TYPE = 2; // euclidean, default
    if (boost::iequals(dist_type, "manhattan")) ANN_DIST_TYPE = 1; // manhattan

    // since KNN search will always return the query point itself, so add 1
    // to make sure returning min_samples number of results
    //min_samples = min_samples + 1;
    GalElement* gal = new GalElement[rows];

    ANNkd_tree* kdTree = new ANNkd_tree(input_data, rows, columns);
    ANNidxArray nnIdx = new ANNidx[knn+1];
    ANNdistArray dists = new ANNdist[knn+1];
    for (int i=0; i<rows; ++i) {
        kdTree->annkSearch(input_data[i], (int)knn+1, nnIdx, dists, eps);
        //core_d[i] = sqrt(dists[min_samples-1]);
        gal[i].SetSizeNbrs(knn);
        for (int j=0; j<knn; j++) {
            gal[i].SetNbr(j, nnIdx[j+1], 1.0);
        }
    }
    delete[] nnIdx;
    delete[] dists;
    delete kdTree;
    for (int i=0; i<rows; i++) delete[] input_data[i];
    delete[] input_data;


    GalWeight* gw = new GalWeight();
    gw->num_obs = rows;
    gw->wflnm = "";
    gw->id_field = "";
    gw->gal = gal;
    gw->GetNbrStats();

    // intersection weights
    std::vector<GeoDaWeight*> two_weights;
    two_weights.push_back(sw);
    two_weights.push_back(gw);
    GalWeight* intersect_w = WeightUtils::WeightsIntersection(two_weights);

    // compute cnbrs (number of common neighbors), p value
    GalElement* new_gal = intersect_w->gal;
    std::vector<double> val_cnbrs(rows);
    for (int i=0; i<rows; ++i) {
        val_cnbrs[i] = (double)new_gal[i].Size();
    }

    // clean up weights
    delete[] new_gal;
    delete gw;
    delete sw;

    int k = (int)knn;
    std::vector<double> pval_dict(knn,  -1);
    for (int v=1; v<k; ++v) {
        // p = C(k,v).C(N-k,k-v) / C(N,k),
        pval_dict[v] = Gda::combinatorial(k, v) * Gda::combinatorial(rows-k-1, k-v);
        pval_dict[v] /= Gda::combinatorial(rows-1, k);
    }
    std::vector<double> val_p(rows);
    for (int i=0; i<rows; ++i) {
        val_p[i] = pval_dict[val_cnbrs[i]];
    }

    // save the results: cnbrs (number of common neighbors), p value
    std::vector<std::vector<double> > result;
    result.push_back(val_cnbrs);
    result.push_back(val_p);

    return result;
}
