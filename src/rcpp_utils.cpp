// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 12/23/2020 init rcpp_utils.cpp

#include <Rcpp.h>
using namespace Rcpp;

#include "libgeoda_src/gda_data.h"

//  [[Rcpp::export]]
Rcpp::List p_gda_demean(Rcpp::List data)
{
  std::vector<std::vector<double> > _data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    _data.push_back(vals);
  }

  std::vector<std::vector<double> > demean_data = gda_demean(_data);

  Rcpp::List out(data.size());
  for (int i=0; i< data.size(); ++i) {
    std::vector<double>& vals = demean_data[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_gda_standardize(Rcpp::List data)
{
  std::vector<std::vector<double> > _data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    _data.push_back(vals);
  }

  std::vector<std::vector<double> > std_data = gda_standardize(_data);

  Rcpp::List out(data.size());
  for (int i=0; i< data.size(); ++i) {
    std::vector<double>& vals = std_data[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
}

//  [[Rcpp::export]]
Rcpp::List p_gda_mad(Rcpp::List data)
{
  std::vector<std::vector<double> > _data;
  for (int i=0; i< data.size(); ++i) {
    Rcpp::NumericVector tmp = data[i];
    std::vector<double> vals = as<std::vector<double> >(tmp);
    _data.push_back(vals);
  }

  std::vector<std::vector<double> > std_data = gda_standardize_mad(_data);

  Rcpp::List out(data.size());
  for (int i=0; i< data.size(); ++i) {
    std::vector<double>& vals = std_data[i];
    Rcpp::NumericVector tmp_vals(vals.begin(), vals.end());
    out[i] = tmp_vals;
  }

  return out;
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
