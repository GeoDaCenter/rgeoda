// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 12/23/2020 init rcpp_lisa.cpp

#include <Rcpp.h>
#include "libgeoda_src/weights/GeodaWeight.h"
#include "libgeoda_src/sa/LISA.h"
#include "libgeoda_src/gda_sa.h"

using namespace Rcpp;

//  [[Rcpp::export]]
void p_LISA__Run(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  // invoke the function
  ptr->Run();
}

//  [[Rcpp::export]]
void p_LISA__SetNumPermutations(SEXP xp, int num_perm)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  // invoke the function
  ptr->SetNumPermutations(num_perm);
}

//  [[Rcpp::export]]
void p_LISA__SetNumThreads(SEXP xp, int num_threads)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  // invoke the function
  ptr->SetNumThreads(num_threads);
}

//  [[Rcpp::export]]
std::vector<double> p_LISA__GetLISAValues(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetLISAValues();
}

//  [[Rcpp::export]]
std::vector<double> p_LISA__GetLocalSignificanceValues(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetLocalSignificanceValues();
}

//  [[Rcpp::export]]
std::vector<int> p_LISA__GetClusterIndicators(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetClusterIndicators();
}

//  [[Rcpp::export]]
std::vector<int> p_LISA__GetNumNeighbors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetNumNeighbors();
}

//  [[Rcpp::export]]
void p_LISA__SetSignificanceCutoff(SEXP xp, double cutoff)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  ptr->SetSignificanceCutoff(cutoff);
}

//  [[Rcpp::export]]
std::vector<std::string> p_LISA__GetLabels(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetLabels();
}

//  [[Rcpp::export]]
std::vector<std::string> p_LISA__GetColors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetColors();
}

//  [[Rcpp::export]]
double p_LISA__GetBO(SEXP xp, double pval)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetBO(pval);
}

//  [[Rcpp::export]]
double p_LISA__GetFDR(SEXP xp, double pval)
{
  // grab the object as a XPtr (smart pointer) to LISA 
  Rcpp::XPtr<LISA> ptr(xp);

  return ptr->GetFDR(pval);
}

//  [[Rcpp::export]]
SEXP p_localmoran(SEXP xp_w, NumericVector data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_localmoran(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_localgeary(SEXP xp_w, NumericVector data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_localgeary(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_localmultigeary(SEXP xp_w, Rcpp::List& data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n_obs = w->GetNumObs();
  int n_vars = data.size();
  std::vector<std::vector<bool> > undefs(n_vars);
  std::vector<std::vector<double> > raw_data(n_vars);

  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i].resize(n_obs);
    undefs[i].resize(n_obs, false);
    for (int j=0; j< n_obs; ++j) {
      raw_data[i][j] = tmp[j];
      undefs[i][j] = undefs[i][j] || tmp.is_na(i);
    }
  }

  LISA* lisa = gda_localmultigeary(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}


//  [[Rcpp::export]]
SEXP p_localg(SEXP xp_w, NumericVector data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_localg(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_localgstar(SEXP xp_w, NumericVector data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_localgstar(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_localjoincount(SEXP xp_w, NumericVector data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_localjoincount(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_localmultijoincount(SEXP xp_w, Rcpp::List& data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n_obs = w->GetNumObs();
  int n_vars = data.size();
  std::vector<std::vector<bool> > undefs(n_vars);
  std::vector<std::vector<double> > raw_data(n_vars);

  for (int i=0; i< n_vars; ++i) {
    Rcpp::NumericVector tmp = data[i];
    raw_data[i].resize(n_obs);
    undefs[i].resize(n_obs, false);
    for (int j=0; j< n_obs; ++j) {
      raw_data[i][j] = tmp[j];
      undefs[i][j] = undefs[i][j] || tmp.is_na(i);
    }
  }

  LISA* lisa = gda_localmultijoincount(w, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_quantilelisa(SEXP xp_w, int k, int quantile, NumericVector& data, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  int n = data.size();
  std::vector<double> raw_data(n);
  std::vector<bool> undefs(n, false);

  for (int i=0; i< data.size(); ++i) {
    raw_data[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  LISA* lisa = gda_quantilelisa(w, k, quantile, raw_data, undefs, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}

//  [[Rcpp::export]]
SEXP p_multiquantilelisa(SEXP xp_w, NumericVector& k_s, NumericVector& q_s, Rcpp::List& data_s, int permutations, double significance_cutoff, int cpu_threads, int seed)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp_w);
  GeoDaWeight* w = static_cast<GeoDaWeight*> (R_ExternalPtrAddr(ptr));

  std::vector<int> ks = as<std::vector<int> >(k_s);
  std::vector<int> qs = as<std::vector<int> >(q_s);

  int n = data_s.size();
  std::vector<std::vector<double> > raw_data_s(n);
  std::vector<std::vector<bool> > undefs_s(n);

  for (int i=0; i< n; ++i) {
    Rcpp::NumericVector tmp = data_s[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);

    raw_data_s[i] = vals;

    for (int j=0; j< tmp.size(); ++j) {
      undefs_s[i].push_back(tmp.is_na(j));
    }
  }

  LISA* lisa = gda_multiquantilelisa(w, ks, qs, raw_data_s, undefs_s, significance_cutoff, cpu_threads, permutations, seed);

  Rcpp::XPtr<LISA> lisa_ptr(lisa, true);
  return lisa_ptr;
}
