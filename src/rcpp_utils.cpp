// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 12/23/2020 init rcpp_utils.cpp

#include <Rcpp.h>
using namespace Rcpp;

#include "libgeoda/gda_data.h"

//  [[Rcpp::export]]
bool p_gda_isbinary(Rcpp::NumericVector& values)
{
  int n = values.size();
  for (int i=0; i<n; ++i) {
    if (values[i] != 0 && values[i] != 1) {
      return false;
    }
  }
  return true;
}



// [[Rcpp::export]]
Rcpp::NumericVector p_naturalbreaks(int k, Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_naturalbreaks(k, vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_quantilebreaks(int k, Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_quantilebreaks(k, vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_hinge15breaks(Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_hinge15breaks(vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_hinge30breaks(Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_hinge30breaks(vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_percentilebreaks(Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_percentilebreaks(vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_stddevbreaks(Rcpp::NumericVector data)
{
  int num_obs = data.size();

  std::vector<double> vals(num_obs);
  std::vector<bool> undefs(num_obs, false);

  for (int i=0; i< num_obs; ++i) {
    vals[i] = data[i];
    undefs[i] = data.is_na(i);
  }

  std::vector<double> result = gda_stddevbreaks(vals, undefs);

  Rcpp::NumericVector out(result.begin(), result.end());
  return out;
}
