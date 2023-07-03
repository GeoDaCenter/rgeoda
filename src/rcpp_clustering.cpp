// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 1/4/2021 init rcpp_clustering.cpp
// 2/11/2021 add scale_method

#include <Rcpp.h>
using namespace Rcpp;

#include "libgeoda/gda_clustering.h"
#include "libgeoda/GenUtils.h"
#include "libgeoda/gda_interface.h"
#include "libgeoda/libgeoda.h"

Rcpp::List _create_clustering_result(int num_obs, const std::vector<std::vector<int> >& cluster_ids,
                                     const std::vector<std::vector<double> >& raw_data)
{
  std::vector<int> clusters = GenUtils::flat_2dclusters(num_obs, cluster_ids);

  double between_ss = gda_betweensumofsquare(cluster_ids, raw_data);
  double total_ss = gda_totalsumofsquare(raw_data);
  double ratio = between_ss / total_ss;
  std::vector<double> within_ss = gda_withinsumofsquare(cluster_ids, raw_data);

  Rcpp::IntegerVector out_clusters(clusters.begin(), clusters.end());
  Rcpp::NumericVector out_withinss(within_ss.begin(), within_ss.end());

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("Clusters") = out_clusters,
    Rcpp::Named("Total sum of squares") = total_ss,
    Rcpp::Named("Within-cluster sum of squares") = out_withinss,
    Rcpp::Named("Total within-cluster sum of squares") = between_ss,
    Rcpp::Named("The ratio of between to total sum of squares") = ratio
  );

  return out;
}

double** rdist_matrix(int n, NumericVector& rdist)
{
  // rdist is stored as lower part triangle, column wise
  if (rdist.size() == 0) return NULL;

  std::vector<double> dist = as<std::vector<double> >(rdist);
  double** matrix = (double**)malloc(n*sizeof(double*));
  matrix[0] = NULL;
  for (int i = 1; i < n; i++) {
    matrix[i] = (double*)malloc(i*sizeof(double));
  }
  int m = (n - 1) * n / 2;
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      int r = i > j ? i : j;
      int c = i < j ? i : j;
      int idx = m - (n - c - 1) * (n - c) / 2 + (r -c) -1 ;
      matrix[i][j] = dist[idx];
    }
  }
  return matrix;
}

//  [[Rcpp::export]]
Rcpp::List p_skater(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data(n_vars);

  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);

  int num_obs = w->GetNumObs();
  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_skater(k, w, raw_data, scale_method, distance_method, raw_bound, min_bound, seed, cpu_threads, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_redcap(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string redcap_method, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data(n_vars);

  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);

  int num_obs = w->GetNumObs();
  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_redcap(k, w, raw_data, scale_method, redcap_method, distance_method, raw_bound, min_bound, seed, cpu_threads, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_schc(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string linkage_method, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data(n_vars);

  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  int num_obs = w->GetNumObs();
  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<std::vector<int> > cluster_ids = gda_schc(k, w, raw_data, scale_method, linkage_method, distance_method, raw_bound, min_bound, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_greedy(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                  int iterations, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_maxp_greedy(w, raw_data, scale_method, iterations, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_sa(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                     int iterations, double cooling_rate, int sa_maxit, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_maxp_sa(w, raw_data, scale_method, iterations, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_tabu(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                       int iterations, int tabu_length, int conv_tabu, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_maxp_tabu(w, raw_data, scale_method, iterations, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_greedy(int p, SEXP xp_w, Rcpp::List& data, int n_vars,  NumericVector& bound_vals, double min_bound, int inits,
                        NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to Weight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  std::vector<std::vector<int> > cluster_ids = gda_azp_greedy(p, w, raw_data, scale_method, inits, min_bounds, max_bounds, raw_init_regions, distance_method, seed, dist_matrix);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_sa(int p, SEXP xp_w, Rcpp::List& data, int n_vars, double cooling_rate, int sa_maxit,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_azp_sa(p, w, raw_data, scale_method, inits, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_tabu(int p, SEXP xp_w, Rcpp::List& data, int n_vars, int tabu_length, int conv_tabu,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, NumericVector& rdist)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = w->GetNumObs();
  std::vector<std::vector<double> > raw_data(n_vars);
  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i] = as<std::vector<double> >(tmp);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions;// = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  double** dist_matrix = rdist_matrix(num_obs, rdist);

  std::vector<std::vector<int> > cluster_ids = gda_azp_tabu(p, w, raw_data, scale_method, inits, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed, dist_matrix);

  if (dist_matrix) {
    for (int i = 1; i < num_obs; i++) {
      free(dist_matrix[i]);
    }
  }
  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_spatialvalidation(SEXP xp_geoda, NumericVector& clusters, SEXP xp_w)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  Rcpp::XPtr<GeoDaWeight> ptr_w(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr_w));

  int n = clusters.size();
  std::vector<int> raw_clusters(n);

  for (int i=0; i< n; ++i) {
    raw_clusters[i] = clusters[i];
  }

  ValidationResult result = gda_spatialvalidation(geoda, raw_clusters, w);

  Rcpp::NumericVector jc_v1, jc_v2, jc_v3, jc_v4, jc_v5;
  int n_clusters = result.joincount_ratio.size();
  for (int i = 0; i < n_clusters; ++i) {
    jc_v1.push_back(i + 1);
    jc_v2.push_back(result.joincount_ratio[i].n);
    jc_v3.push_back(result.joincount_ratio[i].totalNeighbors);
    jc_v4.push_back(result.joincount_ratio[i].totalJoinCount);
    jc_v5.push_back(result.joincount_ratio[i].ratio);
  }
  Rcpp::DataFrame out_joincount = Rcpp::DataFrame::create(
    Rcpp::Named("Cluster") = jc_v1,
    Rcpp::Named("N") = jc_v2,
    Rcpp::Named("Neighbors") = jc_v3,
    Rcpp::Named("Join Count") = jc_v4,
    Rcpp::Named("Ratio") = jc_v5
  );

  JoinCountRatio all_jcr = gda_all_joincount_ratio(result.joincount_ratio);
  Rcpp::List out_all_joincount = Rcpp::List::create(
    Rcpp::Named("N") = all_jcr.n,
    Rcpp::Named("Neighbors") = all_jcr.totalNeighbors,
    Rcpp::Named("Join Count") = all_jcr.totalJoinCount,
    Rcpp::Named("Ratio") = all_jcr.ratio
  );

  Rcpp::List out_fragmentation = Rcpp::List::create(
    Rcpp::Named("#Clusters") = result.fragmentation.n,
    Rcpp::Named("Entropy") = result.fragmentation.entropy,
    Rcpp::Named("Entropy*") = result.fragmentation.std_entropy,
    Rcpp::Named("Simpson") = result.fragmentation.simpson,
    Rcpp::Named("Simpson*") = result.fragmentation.std_simpson
  );

  if (result.spatially_constrained) {
    Rcpp::NumericVector compact_v1, compact_v2, compact_v3, compact_v4;
    n_clusters = result.cluster_compactness.size();
    for (int i = 0; i < n_clusters; ++i) {
      compact_v1.push_back(i + 1);
      compact_v2.push_back(result.cluster_compactness[i].area);
      compact_v3.push_back(result.cluster_compactness[i].perimeter);
      compact_v4.push_back(result.cluster_compactness[i].isoperimeter_quotient);
    }
    Rcpp::DataFrame out_compact = Rcpp::DataFrame::create(
      Rcpp::Named("Cluster") = compact_v1,
      Rcpp::Named("Area") = compact_v2,
      Rcpp::Named("Perimeter") = compact_v3,
      Rcpp::Named("IPC") = compact_v4
    );

    Rcpp::NumericVector diameter_v1, diameter_v2, diameter_v3;
    n_clusters = result.cluster_diameter.size();
    for (int i = 0; i < n_clusters; ++i) {
      diameter_v1.push_back(i + 1);
      diameter_v2.push_back(result.cluster_diameter[i].steps);
      diameter_v3.push_back(result.cluster_diameter[i].ratio);
    }
    Rcpp::DataFrame out_diameter = Rcpp::DataFrame::create(
      Rcpp::Named("Cluster") = diameter_v1,
      Rcpp::Named("Steps") = diameter_v2,
      Rcpp::Named("Ratio") = diameter_v3
    );
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("IsSpatiallyConstrained") = result.spatially_constrained,
      Rcpp::Named("Fragmentation") = out_fragmentation,
      Rcpp::Named("SubclusterFragmentation") = "N/A: clusters are spatially constrained.",
      Rcpp::Named("JoinCountRatio") = out_joincount,
      Rcpp::Named("AllJoinCountRatio") = out_all_joincount,
      Rcpp::Named("Compactness") = out_compact,
      Rcpp::Named("Diameter") = out_diameter
    );

    return out;

  } else {
    Rcpp::NumericVector subfrag_v1, subfrag_v2, subfrag_v3, subfrag_v4, subfrag_v5;
    Rcpp::NumericVector subfrag_v6, subfrag_v7, subfrag_v8, subfrag_v9, subfrag_v10, subfrag_v11;
    n_clusters = result.cluster_fragmentation.size();
    for (int i = 0; i < n_clusters; ++i) {
      subfrag_v1.push_back(i + 1);
      subfrag_v2.push_back(result.joincount_ratio[i].n);
      subfrag_v3.push_back(result.cluster_fragmentation[i].fraction);
      subfrag_v4.push_back(result.cluster_fragmentation[i].n);
      subfrag_v5.push_back(result.cluster_fragmentation[i].entropy);
      subfrag_v6.push_back(result.cluster_fragmentation[i].std_entropy);
      subfrag_v7.push_back(result.cluster_fragmentation[i].simpson);
      subfrag_v8.push_back(result.cluster_fragmentation[i].std_simpson);
      subfrag_v9.push_back(result.cluster_fragmentation[i].min_cluster_size);
      subfrag_v10.push_back(result.cluster_fragmentation[i].max_cluster_size);
      subfrag_v11.push_back(result.cluster_fragmentation[i].mean_cluster_size);
    }
    Rcpp::DataFrame out_subfrag = Rcpp::DataFrame::create(
      Rcpp::Named("Cluster") = subfrag_v1,
      Rcpp::Named("N") = subfrag_v2,
      Rcpp::Named("Fraction") = subfrag_v3,
      Rcpp::Named("#Sub") = subfrag_v4,
      Rcpp::Named("Entropy") = subfrag_v5,
      Rcpp::Named("Entropy*") = subfrag_v6,
      Rcpp::Named("Simpson") = subfrag_v7,
      Rcpp::Named("Simpson*") = subfrag_v8,
      Rcpp::Named("Min") = subfrag_v9,
      Rcpp::Named("Max") = subfrag_v10,
      Rcpp::Named("Mean") = subfrag_v11
    );

    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("IsSpatiallyConstrained") = result.spatially_constrained,
      Rcpp::Named("Fragmentation") = out_fragmentation,
      Rcpp::Named("SubclusterFragmentation") = out_subfrag,
      Rcpp::Named("JoinCountRatio") = out_joincount,
      Rcpp::Named("AllJoinCountRatio") = out_all_joincount,
      Rcpp::Named("Compactness") = "N/A: clusters are not spatially constrained.",
      Rcpp::Named("Diameter") = "N/A: clusters are not spatially constrained."
    );

    return out;
  }

}

//  [[Rcpp::export]]
Rcpp::List p_joincount_ratio(NumericVector& clusters, SEXP xp_w)
{
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<GeoDaWeight> ptr_w(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr_w));

  int n = clusters.size();
  std::vector<int> raw_clusters(n);

  for (int i=0; i< n; ++i) {
    raw_clusters[i] = clusters[i];
  }

  std::vector<JoinCountRatio> items = gda_joincount_ratio(raw_clusters, w);
  JoinCountRatio all_jcr = gda_all_joincount_ratio(items);

  Rcpp::NumericVector jc_v1, jc_v2, jc_v3, jc_v4, jc_v5;
  int n_clusters = items.size();
  for (int i = 0; i < n_clusters; ++i) {
    jc_v1.push_back(i + 1);
    jc_v2.push_back(items[i].n);
    jc_v3.push_back(items[i].totalNeighbors);
    jc_v4.push_back(items[i].totalJoinCount);
    jc_v5.push_back(items[i].ratio);
  }
  Rcpp::DataFrame out_joincount = Rcpp::DataFrame::create(
    Rcpp::Named("Cluster") = jc_v1,
    Rcpp::Named("N") = jc_v2,
    Rcpp::Named("Neighbors") = jc_v3,
    Rcpp::Named("Join Count") = jc_v4,
    Rcpp::Named("Ratio") = jc_v5
  );

  Rcpp::List out_all_joincount = Rcpp::List::create(
    Rcpp::Named("N") = all_jcr.n,
    Rcpp::Named("Neighbors") = all_jcr.totalNeighbors,
    Rcpp::Named("Join Count") = all_jcr.totalJoinCount,
    Rcpp::Named("Ratio") = all_jcr.ratio
  );

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("JoinCountRatio") = out_joincount,
    Rcpp::Named("AllJoinCountRatio") = out_all_joincount
  );

  return out;
}

//  [[Rcpp::export]]
Rcpp::NumericVector p_make_spatial(NumericVector& clusters, SEXP xp_w)
{
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<GeoDaWeight> ptr_w(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr_w));

  int n = clusters.size();
  std::vector<int> raw_clusters(n);

  for (int i=0; i< n; ++i) {
    raw_clusters[i] = clusters[i];
  }

  std::vector<int> result = gda_makespatial(raw_clusters, w);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}
