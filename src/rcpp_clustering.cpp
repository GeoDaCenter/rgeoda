// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 1/4/2021 init rcpp_clustering.cpp

#include <Rcpp.h>
using namespace Rcpp;

#include "libgeoda_src/gda_clustering.h"

//  [[Rcpp::export]]
double p_betweensumofsquare(Rcpp::List& solution, Rcpp::List& data)
{
  std::vector<std::vector<int> > raw_sol;
  std::vector<std::vector<double> > raw_data;

  for (int i=0; i< solution.size(); ++i) {
    Rcpp::NumericVector tmp = solution[i];
    std::vector<int> vals = as<std::vector<int> >(tmp);
    raw_sol.push_back(vals);
  }

  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  double out = gda_betweensumofsquare(raw_sol, raw_data);
  return out;
}

//  [[Rcpp::export]]
double p_totalsumofsquare(Rcpp::List& data)
{
  std::vector<std::vector<double> > raw_data;

  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  double out = gda_totalsumofsquare(raw_data);
  return out;
}

//  [[Rcpp::export]]
double p_withinsumofsquare(Rcpp::List& solution, Rcpp::List& data)
{
  std::vector<std::vector<int> > raw_sol;
  std::vector<std::vector<double> > raw_data;

  for (int i=0; i< solution.size(); ++i) {
    Rcpp::NumericVector tmp = solution[i];
    std::vector<int> vals = as<std::vector<int> >(tmp);
    raw_sol.push_back(vals);
  }

  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  double out = gda_withinsumofsquare(raw_sol, raw_data);
  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_skater(int k, SEXP xp_w, Rcpp::List& data, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<std::vector<int> > clusters = gda_skater(k, w, raw_data, distance_method, raw_bound, min_bound, seed, cpu_threads);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_redcap(int k, SEXP xp_w, Rcpp::List& data, std::string redcap_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound, int seed, int cpu_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<std::vector<int> > clusters = gda_redcap(k, w, raw_data, redcap_method, distance_method, raw_bound, min_bound, seed, cpu_threads);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_schc(int k, SEXP xp_w, Rcpp::List& data, std::string linkage_method, std::string distance_method,
                    NumericVector& bound_vals, double min_bound)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<std::vector<int> > clusters = gda_schc(k, w, raw_data, linkage_method, distance_method, raw_bound, min_bound);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_greedy(SEXP xp_w, Rcpp::List& data, NumericVector& bound_vals, double min_bound,
                  int iterations, NumericVector& init_regions, std::string distance_method, int seed, int cpu_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_maxp_greedy(w, raw_data, iterations, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_sa(SEXP xp_w, Rcpp::List& data, NumericVector& bound_vals, double min_bound,
                     int iterations, double cooling_rate, int sa_maxit, NumericVector& init_regions, std::string distance_method, int seed, int cpu_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_maxp_sa(w, raw_data, iterations, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_maxp_tabu(SEXP xp_w, Rcpp::List& data, NumericVector& bound_vals, double min_bound,
                       int iterations, int tabu_length, int conv_tabu, NumericVector& init_regions, std::string distance_method, int seed, int cpu_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_maxp_tabu(w, raw_data, iterations, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed, cpu_threads);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_azp_greedy(int p, SEXP xp_w, Rcpp::List& data, NumericVector& bound_vals, double min_bound, int inits,
                        NumericVector& init_regions, std::string distance_method, int seed)
{
  // grab the object as a XPtr (smart pointer) to Weight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_azp_greedy(p, w, raw_data, inits, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_azp_sa(int p, SEXP xp_w, Rcpp::List& data, double cooling_rate, int sa_maxit,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string distance_method, int seed)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_azp_sa(p, w, raw_data, inits, cooling_rate, sa_maxit, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_azp_tabu(int p, SEXP xp_w, Rcpp::List& data, int tabu_length, int conv_tabu,
                    NumericVector& bound_vals, double min_bound, int inits,
                    NumericVector& init_regions, std::string distance_method, int seed)
{
  // grab the object as a XPtr (smart pointer) to LISA
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int num_obs = 0;
  std::vector<std::vector<double> > raw_data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    num_obs = vals.size();
    raw_data.push_back(vals);
  }

  std::vector<double> raw_bound = as<std::vector<double> >(bound_vals);
  std::vector<int> raw_init_regions;// = as<std::vector<int> >(init_regions);

  std::vector<std::pair<double, std::vector<double> > > min_bounds, max_bounds;
  if (raw_bound.size() == num_obs) {
    min_bounds.push_back(std::make_pair(min_bound, raw_bound));
  }

  std::vector<std::vector<int> > clusters = gda_azp_tabu(p, w, raw_data, inits, tabu_length, conv_tabu, min_bounds, max_bounds, raw_init_regions, distance_method, seed);

  Rcpp::List out(clusters.size());
  for (int i=0; i< clusters.size(); ++i) {
    std::vector<int>& vals = clusters[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}
