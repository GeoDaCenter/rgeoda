//
// Created by Xun Li on 1/6/21.
//

#ifndef GEODA_AZP_WRAPPER_H
#define GEODA_AZP_WRAPPER_H


#include <vector>

class GeoDa;
class GalElement;
class GeoDaWeight;
class ZoneControl;
class RegionMaker;
class RawDistMatrix;

class azp_wrapper {
public:
    azp_wrapper(int p, GeoDaWeight *w,
                 const std::vector<std::vector<double> >& data,
                 int inits,
                 const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                 const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                 const std::vector<int>& init_regions,
                 const std::string &distance_method,
                 int rnd_seed);

    virtual ~azp_wrapper();

    virtual const std::vector<std::vector<int> > GetClusters();

    virtual RegionMaker* RunAZP() = 0;

    virtual void CreateController(const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds);

    virtual void Run();

protected:
    int p;

    int num_obs;

    int n_cols;

    int inits;

    std::string distance_method;

    std::vector<std::vector<double> > data;

    GalElement *gal;

    double **input_data;

    RawDistMatrix *dm;

    std::vector<ZoneControl> controllers;

    std::vector<int> init_regions;

    int rnd_seed;

    std::vector<std::vector<int> > cluster_ids;
};

class azp_greedy_wrapper : public azp_wrapper {
public:
    azp_greedy_wrapper(int p, GeoDaWeight *w,
                        const std::vector<std::vector<double> >& data,
                        int inits,
                        const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                        const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                        const std::vector<int>& init_regions,
                        const std::string &distance_method,
                        int rnd_seed);

    virtual ~azp_greedy_wrapper();

    virtual RegionMaker* RunAZP();

protected:
    double cooling_rate;

    int sa_maxit;
};

class azp_sa_wrapper : public azp_wrapper {
public:
    azp_sa_wrapper(int p, GeoDaWeight *w,
                    const std::vector<std::vector<double> >& data,
                    int inits,
                    double cooling_rate,
                    int sa_maxit,
                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                    const std::vector<int>& init_regions,
                    const std::string &distance_method,
                    int rnd_seed);

    virtual ~azp_sa_wrapper();

    virtual RegionMaker* RunAZP();

protected:
    double cooling_rate;

    int sa_maxit;
};

class azp_tabu_wrapper : public azp_wrapper {
public:
    azp_tabu_wrapper(int p, GeoDaWeight *w,
                    const std::vector<std::vector<double> >& data,
                    int inits,
                    int tabu_length,
                    int conv_tabu,
                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                    const std::vector<int>& init_regions,
                    const std::string &distance_method,
                    int rnd_seed);

    virtual ~azp_tabu_wrapper();

    virtual RegionMaker* RunAZP();

protected:
    int tabu_length;

    int conv_tabu;
};
#endif //GEODA_azp_wrapper_H
