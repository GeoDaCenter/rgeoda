//
// Created by Xun Li on 9/26/19.
//

#include <boost/algorithm/string.hpp>
#include "../weights/GalWeight.h"
#include "../GenUtils.h"
#include "cluster.h"
#include "azp.h"
#include "maxp_wrapper.h"

maxp_wrapper::~maxp_wrapper() {

}

maxp_wrapper::maxp_wrapper(GeoDaWeight *w,
                           const std::vector<std::vector<double> >& data,
                           int iterations,
                           const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                           const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                           const std::vector<int>& init_regions,
                           const std::string &distance_method,
                           int rnd_seed,
                           int cpu_threads)
                           : num_obs(w->num_obs), n_cols(data.size()), iterations(iterations),
                           distance_method(distance_method), data(data), init_regions(init_regions), rnd_seed(rnd_seed),
                           cpu_threads(cpu_threads)
{
    gal = Gda::GetGalElement(w);

    // create bounds
    CreateController(min_bounds, max_bounds);
}

void maxp_wrapper::Run() {
        if (gal) {

            // get distance matrix
            input_data = new double*[num_obs];
            int** mask = new int*[num_obs];
            for (int i=0; i<num_obs; ++i) {
                input_data[i] = new double[n_cols];
                mask[i] = new int[n_cols];
                for (int j=0; j<n_cols; ++j) mask[i][j] = 1.0;
            }
            for (int i=0; i<n_cols; ++i) {
                // the data will be standardized in the caller
                //std::vector<double>& vals = data[i];
                //GenUtils::StandardizeData(vals);
                for (int r=0; r<num_obs; ++r) {
                    input_data[r][i] = data[i][r];
                }
            }

            char dist = 'e';
            if (boost::iequals(distance_method, "manhattan")) dist = 'b';
            int transpose = 0; // row wise
            double* weight = new double[n_cols];
            for (int i=0; i<n_cols; ++i) weight[i] = 1.0;

            double** ragged_distances = distancematrix(num_obs, n_cols, input_data,  mask, weight, dist, transpose);
            dm = new RawDistMatrix(ragged_distances);

            MaxpRegion* maxp = RunMaxp();

            std::vector<int> final_solution = maxp->GetResults();

            delete maxp;

            std::map<int, std::vector<int> > solution;
            for (size_t i=0; i<final_solution.size(); ++i) {
                solution[final_solution[i]].push_back(i);
            }
            std::map<int, std::vector<int> >::iterator it;
            for (it = solution.begin(); it != solution.end(); ++it) {
                cluster_ids.push_back(it->second);
            }

            for (int i = 1; i < num_obs; i++) free(ragged_distances[i]);
            free(ragged_distances);

            delete dm;
        }
}

void maxp_wrapper::CreateController(const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds)
{
    // min bounds
    for (size_t i=0; i<min_bounds.size(); ++i) {
        const std::pair<double, std::vector<double> >& bound = min_bounds[i];
        double min_bound = bound.first;
        std::vector<double> bound_vals = bound.second;
        ZoneControl zc(bound_vals);
        zc.AddControl(ZoneControl::SUM,
                      ZoneControl::MORE_THAN, min_bound);
        controllers.push_back(zc);
    }
    // max bounds
    for (size_t i=0; i<max_bounds.size(); ++i) {
        const std::pair<double, std::vector<double> >& bound = max_bounds[i];
        double max_bound = bound.first;
        std::vector<double> bound_vals = bound.second;
        ZoneControl zc(bound_vals);
        zc.AddControl(ZoneControl::SUM,
                      ZoneControl::LESS_THAN, max_bound);
        controllers.push_back(zc);
    }
}

const std::vector<std::vector<int> > maxp_wrapper::GetClusters() {
    return cluster_ids;
}


maxp_greedy_wrapper::~maxp_greedy_wrapper(){

}

maxp_greedy_wrapper::maxp_greedy_wrapper(GeoDaWeight *w,
                                         const std::vector<std::vector<double> >& data,
                                         int iterations,
                                         const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                         const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                         const std::vector<int>& init_regions,
                                         const std::string &distance_method,
                                         int rnd_seed,
                                         int cpu_threads)
        : maxp_wrapper(w, data, iterations, min_bounds, max_bounds, init_regions, distance_method, rnd_seed, cpu_threads)
{
    Run();
}

MaxpRegion* maxp_greedy_wrapper::RunMaxp() {
    MaxpRegion* maxp = new MaxpGreedy(iterations, gal, input_data, dm, num_obs, n_cols, controllers,
            init_regions,rnd_seed, cpu_threads);

    return maxp;
}

maxp_sa_wrapper::~maxp_sa_wrapper(){

}

maxp_sa_wrapper::maxp_sa_wrapper(GeoDaWeight *w,
                                 const std::vector<std::vector<double> >& data,
                                 int iterations,
                                 double cooling_rate,
                                 int sa_maxit,
                                 const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                 const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                 const std::vector<int>& init_regions,
                                 const std::string &distance_method,
                                 int rnd_seed,
                                 int cpu_threads)
: maxp_wrapper(w, data, iterations, min_bounds, max_bounds, init_regions, distance_method, rnd_seed, cpu_threads),
cooling_rate(cooling_rate), sa_maxit(sa_maxit)
{
    Run();
}

MaxpRegion* maxp_sa_wrapper::RunMaxp() {
    MaxpRegion* maxp = new MaxpSA(iterations, gal, input_data, dm, num_obs, n_cols, controllers, cooling_rate, sa_maxit,
            init_regions, rnd_seed, cpu_threads);

    return maxp;
}

maxp_tabu_wrapper::~maxp_tabu_wrapper(){

}

maxp_tabu_wrapper::maxp_tabu_wrapper(GeoDaWeight *w,
                                     const std::vector<std::vector<double> >& data,
                                     int iterations,
                                     int tabu_length,
                                     int conv_tabu,
                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                     const std::vector<int>& init_regions,
                                     const std::string &distance_method,
                                     int rnd_seed,
                                     int cpu_threads)
        : maxp_wrapper(w, data, iterations, min_bounds, max_bounds, init_regions, distance_method, rnd_seed, cpu_threads),
          tabu_length(tabu_length), conv_tabu(conv_tabu)
{
    Run();
}

MaxpRegion* maxp_tabu_wrapper::RunMaxp() {
    MaxpRegion* maxp = new MaxpTabu(iterations, gal, input_data, dm, num_obs, n_cols, controllers, tabu_length,
            conv_tabu, init_regions, rnd_seed, cpu_threads);

    return maxp;
}
