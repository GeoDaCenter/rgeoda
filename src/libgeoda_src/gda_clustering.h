#ifndef __JSGEODSA_GDA_CLUSTERING__
#define __JSGEODSA_GDA_CLUSTERING__

#include <vector>
#include <string>

#include "./weights/GeodaWeight.h"

// APIs of clustering
/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
const std::vector<std::vector<int> > gda_azp_greedy(int p, GeoDaWeight *w,
                                                    const std::vector<std::vector<double> > &data,
                                                    int inits,
                                                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                    const std::vector<int>& init_regions,
                                                    const std::string &distance_method,
                                                    int rnd_seed);

/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param cooling_rate
 * @param sa_maxit
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
const std::vector<std::vector<int> > gda_azp_sa(int p, GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                int inits,
                                                double cooling_rate,
                                                int sa_maxit,
                                                const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                const std::vector<int>& init_regions,
                                                const std::string &distance_method,
                                                int rnd_seed);

/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param tabu_length
 * @param conv_tabu
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
const std::vector<std::vector<int> > gda_azp_tabu(int p, GeoDaWeight *w,
                                                  const std::vector<std::vector<double> > &data,
                                                  int inits,
                                                  int tabu_length,
                                                  int conv_tabu,
                                                  const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                  const std::vector<int>& init_regions,
                                                  const std::string &distance_method,
                                                  int rnd_seed);

/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_maxp_greedy(GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &data,
                                                     int iterations,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed,
                                                     int cpu_threads);


/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param cooling_rate
 * @param sa_maxit
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
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
                                                 int cpu_threads);

/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param tabu_length
 * @param conv_tabu
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
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
                                                   int cpu_threads);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param redcap_method
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_redcap(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &redcap_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_skater(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads);


/**
 *
 * @param k
 * @param w
 * @param data
 * @param linkage_method
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_schc(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &linkage_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound);

/**
 *
 * @param vals
 * @return
 */
double gda_sumofsquares(const std::vector<double>& vals);

/**
 *
 * @param vals
 * @return
 */
double gda_totalsumofsquare(const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param vals
 * @return
 */
double gda_withinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param data
 * @return
 */
double gda_betweensumofsquare(const std::vector<std::vector<int> >& solution,
                              const std::vector<std::vector<double> >& data);

#endif

