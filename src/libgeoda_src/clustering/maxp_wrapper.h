//
// Created by Xun Li on 9/26/19.
//

#ifndef GEODA_MAXP_WRAPPER_H
#define GEODA_MAXP_WRAPPER_H


#include <vector>

class GeoDa;
class GalElement;
class GeoDaWeight;
class ZoneControl;
class MaxpRegion;
class RawDistMatrix;

class maxp_wrapper {
public:
    maxp_wrapper(GeoDaWeight *w,
                 const std::vector<std::vector<double> >& data,
                 int iterations,
                 const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                 const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                 const std::vector<int>& init_regions,
                 const std::string &distance_method,
                 int rnd_seed,
                 int cpu_threads);

    virtual ~maxp_wrapper();

    virtual const std::vector<std::vector<int> > GetClusters();

    virtual MaxpRegion* RunMaxp() = 0;

    virtual void CreateController(const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds);

    virtual void Run();

protected:
    int num_obs;

    int n_cols;

    int iterations;

    std::string distance_method;

    std::vector<std::vector<double> > data;

    GalElement *gal;

    double **input_data;

    RawDistMatrix *dm;

    std::vector<ZoneControl> controllers;

    std::vector<int> init_regions;

    int rnd_seed;

    std::vector<std::vector<int> > cluster_ids;

    int cpu_threads;
};

class maxp_greedy_wrapper : public maxp_wrapper {
public:
    maxp_greedy_wrapper(GeoDaWeight *w,
                        const std::vector<std::vector<double> >& data,
                        int iterations,
                        const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                        const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                        const std::vector<int>& init_regions,
                        const std::string &distance_method,
                        int rnd_seed,
                        int cpu_threads);

    virtual ~maxp_greedy_wrapper();

    virtual MaxpRegion* RunMaxp();

protected:
    double cooling_rate;

    int sa_maxit;
};

class maxp_sa_wrapper : public maxp_wrapper {
public:
    maxp_sa_wrapper(GeoDaWeight *w,
                    const std::vector<std::vector<double> >& data,
                    int iterations,
                    double cooling_rate,
                    int sa_maxit,
                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                    const std::vector<int>& init_regions,
                    const std::string &distance_method,
                    int rnd_seed,
                    int cpu_threads);

    virtual ~maxp_sa_wrapper();

    virtual MaxpRegion* RunMaxp();

protected:
    double cooling_rate;

    int sa_maxit;
};

class maxp_tabu_wrapper : public maxp_wrapper {
public:
    maxp_tabu_wrapper(GeoDaWeight *w,
                    const std::vector<std::vector<double> >& data,
                    int iterations,
                    int tabu_length,
                    int conv_tabu,
                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                    const std::vector<int>& init_regions,
                    const std::string &distance_method,
                    int rnd_seed,
                    int cpu_threads);

    virtual ~maxp_tabu_wrapper();

    virtual MaxpRegion* RunMaxp();

protected:
    int tabu_length;

    int conv_tabu;
};
#endif //GEODA_MAXP_WRAPPER_H
