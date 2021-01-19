#include <iostream>
#include <boost/algorithm/string.hpp>
#include "weights/GeodaWeight.h"
#include "clustering/maxp_wrapper.h"
#include "clustering/redcap_wrapper.h"
#include "clustering/azp_wrapper.h"
#include "clustering/schc_wrapper.h"
#include "GenUtils.h"
#include "gda_clustering.h"


const std::vector<std::vector<int> > gda_azp_greedy(int p, GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &data,
                                                     int inits,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    azp_greedy_wrapper azp(p, w, data, inits, min_bounds, max_bounds, init_regions, distance_method,
            rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_azp_sa(int p, GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                int inits,
                                                double cooling_rate,
                                                int sa_maxit,
                                                const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                const std::vector<int>& init_regions,
                                                const std::string &distance_method,
                                                int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    azp_sa_wrapper azp(p, w, data, inits, cooling_rate, sa_maxit, min_bounds, max_bounds, init_regions, distance_method,
                           rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_azp_tabu(int p, GeoDaWeight *w,
                                                  const std::vector<std::vector<double> > &data,
                                                  int inits,
                                                  int tabu_length,
                                                  int conv_tabu,
                                                  const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                  const std::vector<int>& init_regions,
                                                  const std::string &distance_method,
                                                  int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    azp_tabu_wrapper azp(p, w, data, inits, tabu_length, conv_tabu, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_greedy(GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &data,
                                                     int iterations,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed,
                                                     int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    maxp_greedy_wrapper maxp(w, data, iterations, min_bounds, max_bounds, init_regions, distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_sa(GeoDaWeight *w,
                                              const std::vector<std::vector<double> > &data,
                                              int iterations,
                                              double cooling_rate,
                                              int sa_maxit,
                                              const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                              const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                              const std::vector<int>& init_regions,
                                              const std::string &distance_method,
                                              int rnd_seed,
                                              int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    maxp_sa_wrapper maxp(w, data, iterations, cooling_rate, sa_maxit, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_tabu(GeoDaWeight *w,
                                                   const std::vector<std::vector<double> > &data,
                                                   int iterations,
                                                   int tabu_length,
                                                   int conv_tabu,
                                                   const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                   const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                   const std::vector<int>& init_regions,
                                                   const std::string &distance_method,
                                                   int rnd_seed,
                                                   int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    maxp_tabu_wrapper maxp(w, data, iterations, tabu_length, conv_tabu, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_redcap(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &redcap_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads)
{
    std::vector<std::vector<int> > result;
    unsigned int method = 0;
    if (boost::iequals(redcap_method, "firstorder-singlelinkage")) {
        method = 0;
    } else if  (boost::iequals(redcap_method, "fullorder-completelinkage")){
        method = 1;
    } else if  (boost::iequals(redcap_method, "fullorder-averagelinkage")) {
        method = 2;
    } else if  (boost::iequals(redcap_method, "fullorder-singlelinkage")) {
        method = 3;
    } else if  (boost::iequals(redcap_method, "fullorder-wardlinkage")) {
        method = 4;
    }

    if (w == 0 ||  method > 4) return result;

    if ((int)k > w->num_obs) return result;

    redcap_wrapper redcap(k, w, data, method, distance_method, bound_vals, min_bound, rand_seed, cpu_threads);
    return redcap.GetClusters();
}

const std::vector<std::vector<int> > gda_skater(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads)
{
    return gda_redcap(k, w, data, "firstorder-singlelinkage", distance_method, bound_vals, min_bound, rand_seed, cpu_threads);
}

const std::vector<std::vector<int> > gda_schc(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &linkage_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound)
{
    std::vector<std::vector<int> > result;
    unsigned int method = 0;
    if (boost::iequals(linkage_method, "single")) {
        method = 0;
    } else if  (boost::iequals(linkage_method, "complete")){
        method = 1;
    } else if  (boost::iequals(linkage_method, "average")) {
        method = 2;
    } else if  (boost::iequals(linkage_method, "ward")) {
        method = 3;
    }

    if (w == 0 ||  method > 4) return result;

    if ((int)k > w->num_obs) return result;

    schc_wrapper schc(k, w, data, method, distance_method, bound_vals, min_bound);
    return schc.GetClusters();
}

double gda_sumofsquares(const std::vector<double>& vals)
{
    std::vector<double> data = vals;
    return  GenUtils::SumOfSquares(data);
}

double gda_totalsumofsquare(const std::vector<std::vector<double> >& vals)
{
    double ssq = 0;
    for (size_t i=0; i<vals.size(); ++i) {
        std::vector<double> data = vals[i];
        GenUtils::StandardizeData(data);
        double ss = gda_sumofsquares(data);
        ssq += ss;
    }
    return ssq;
}

double gda_withinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& _data)
{
    double ssq = 0;
    size_t cols = _data.size();

    // standardize data
    std::vector<std::vector<double> > data(cols);
    for (size_t c=0; c<cols; ++c) {
        data[c] = _data[c];
        GenUtils::StandardizeData(data[c]);
    }

    for (size_t c=0; c<cols; ++c) {
        for (size_t i=0; i<solution.size(); ++i) {
            std::vector<double> vals;
            for (size_t j = 0; j < solution[i].size(); ++j) {
                size_t r = solution[i][j];
                vals.push_back(data[c][r]);
            }
            double ss = gda_sumofsquares(vals);
            ssq += ss;
        }
    }
    return ssq;
}

double gda_betweensumofsquare(const std::vector<std::vector<int> >& solution,
                              const std::vector<std::vector<double> >& data)
{
    double totss = gda_totalsumofsquare(data);
    double totwithiness = gda_withinsumofsquare(solution, data);
    double betweenss = totss - totwithiness;
    return betweenss;
}
