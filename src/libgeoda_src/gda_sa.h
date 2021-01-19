#ifndef __JSGEODSA_GDA_SA__
#define __JSGEODSA_GDA_SA__

#include <string>
#include <vector>

class GeoDaWeight;
class LISA;
class BatchLISA;

// APIs of local spatial autocorrelation
/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localmoran(GeoDaWeight *w,
                     const std::vector<double> &data,
                     const std::vector<bool> &undefs, double significance_cutoff,
                     int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
BatchLISA *gda_batchlocalmoran(GeoDaWeight *w, const std::vector<std::vector<double> > &data,
                               const std::vector<std::vector<bool> > &undefs, double significance_cutoff,
                               int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localgeary(GeoDaWeight *w, const std::vector<double> &data,
                     const std::vector<bool> &undefs, double significance_cutoff,
                     int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localmultigeary(GeoDaWeight *w, const std::vector<std::vector<double> > &data,
                     const std::vector<std::vector<bool> > &undefs, double significance_cutoff,
                     int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localjoincount(GeoDaWeight *w, const std::vector<double> &data,
                    const std::vector<bool> &undefs, double significance_cutoff,
                    int nCPUs, int permutations, int last_seed_used);

LISA *gda_localmultijoincount(GeoDaWeight *w, const std::vector<std::vector<double> > &data,
                         const std::vector<std::vector<bool> > &undefs, double significance_cutoff,
                         int nCPUs, int permutations, int last_seed_used);
/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localg(GeoDaWeight *w, const std::vector<double> &data,
                 const std::vector<bool> &undefs, double significance_cutoff,
                 int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_localgstar(GeoDaWeight *w, const std::vector<double> &data,
                     const std::vector<bool> &undefs, double significance_cutoff,
                     int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param quantile
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_quantilelisa(GeoDaWeight *w, unsigned int k, unsigned int quantile, const std::vector<double> &data,
                       const std::vector<bool> &undefs, double significance_cutoff,
                       int nCPUs, int permutations, int last_seed_used);

/**
 *
 * @param w
 * @param k
 * @param quantile
 * @param data
 * @param undefs
 * @return
 */
LISA *gda_multiquantilelisa(GeoDaWeight *w, const std::vector<int>& k_s, const std::vector<int>& quantile_s, const std::vector<std::vector<double> > &data_s,
                       const std::vector<std::vector<bool> > &undefs_s, double significance_cutoff,
                       int nCPUs, int permutations, int last_seed_used);
/**
 *
 * @param lisa
 * @param current_p
 * @return
 */
double gda_fdr(LISA *lisa, double current_p);

/**
 *
 * @param lisa
 * @param current_p
 * @return
 */
double gda_bo(LISA *lisa, double current_p);

#endif
