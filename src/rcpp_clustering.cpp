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

//  [[Rcpp::export]]
Rcpp::List p_skater(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads)
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

  std::vector<std::vector<int> > cluster_ids = gda_skater(k, w, raw_data, scale_method, distance_method, raw_bound, min_bound, seed, cpu_threads);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_redcap(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string redcap_method, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads)
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
  std::vector<std::vector<int> > cluster_ids = gda_redcap(k, w, raw_data, scale_method, redcap_method, distance_method, raw_bound, min_bound, seed, cpu_threads);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_schc(int k, SEXP xp_w, Rcpp::List& data, int n_vars, std::string linkage_method, std::string scale_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound)
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
  std::vector<std::vector<int> > cluster_ids = gda_schc(k, w, raw_data, scale_method, linkage_method, distance_method, raw_bound, min_bound);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_greedy(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                  int iterations, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads)
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

  std::vector<std::vector<int> > cluster_ids = gda_maxp_greedy(w, raw_data, scale_method, iterations, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_sa(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                     int iterations, double cooling_rate, int sa_maxit, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads)
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

  std::vector<std::vector<int> > cluster_ids = gda_maxp_sa(w, raw_data, scale_method, iterations, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_tabu(SEXP xp_w, Rcpp::List& data, int n_vars, NumericVector& bound_vals, double min_bound,
                       int iterations, int tabu_length, int conv_tabu, NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed, int cpu_threads)
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

  std::vector<std::vector<int> > cluster_ids = gda_maxp_tabu(w, raw_data, scale_method, iterations, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_greedy(int p, SEXP xp_w, Rcpp::List& data, int n_vars,  NumericVector& bound_vals, double min_bound, int inits,
                        NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed)
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

  std::vector<std::vector<int> > cluster_ids = gda_azp_greedy(p, w, raw_data, scale_method, inits, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_sa(int p, SEXP xp_w, Rcpp::List& data, int n_vars, double cooling_rate, int sa_maxit,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed)
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

  std::vector<std::vector<int> > cluster_ids = gda_azp_sa(p, w, raw_data, scale_method, inits, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}

//  [[Rcpp::export]]
Rcpp::List p_azp_tabu(int p, SEXP xp_w, Rcpp::List& data, int n_vars, int tabu_length, int conv_tabu,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string scale_method, std::string distance_method, int seed)
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

  std::vector<std::vector<int> > cluster_ids = gda_azp_tabu(p, w, raw_data, scale_method, inits, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  return _create_clustering_result(w->GetNumObs(), cluster_ids, raw_data);
}
