// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
// Author: lixun910@gmail.com
// Changes:
// 10/29/2020 init file

#include <Rcpp.h>
#include "libgeoda_src/libgeoda.h"

using namespace Rcpp;

//  [[Rcpp::export]]
SEXP p_GeoDa__new(std::string file_path)
{
  // create a pointer to an GeoDa object and wrap it
  // as an external pointer
  Rcpp::XPtr<GeoDa> ptr( new GeoDa(file_path.c_str()), true );

  // return the external pointer to the R side
  return ptr;
}

//  [[Rcpp::export]]
SEXP p_GeoDa__new1(std::string layer_name, std::string map_type, int num_features, RawVector wkbs, NumericVector wkb_bytes_len)
{
  std::vector<unsigned char> _wkbs = as<std::vector<unsigned char> >(wkbs);
  std::vector<int> _wkb_bytes_len = as<std::vector<int> >(wkb_bytes_len);

  // create a pointer to an GeoDa object and wrap it
  // as an external pointer
  Rcpp::XPtr<GeoDa> ptr( new GeoDa(layer_name, map_type, num_features, _wkbs.data(), _wkb_bytes_len),
                         true );

  // return the external pointer to the R side
  return ptr;
}

//  [[Rcpp::export]]
int p_GeoDa__GetNumObs(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  int num_obs = ptr->GetNumObs();

  return num_obs;
}


//  [[Rcpp::export]]
int p_GeoDa__GetNumCols(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  int num_cols = ptr->GetNumCols();

  return num_cols;
}

//  [[Rcpp::export]]
int p_GeoDa__GetMapType(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  int map_type = ptr->GetMapType();

  return map_type;
}

//  [[Rcpp::export]]
StringVector p_GeoDa__GetFieldNames(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<std::string> field_nms = ptr->GetFieldNames();

  int n_fields = field_nms.size();

  // convert to Rcpp::StringVector
  StringVector sv_field_nms(n_fields);

  for (int i=0; i<n_fields; ++i) {
    sv_field_nms[i] = field_nms[i];
  }

  return sv_field_nms;
}

//  [[Rcpp::export]]
StringVector p_GeoDa__GetFieldTypes(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<std::string> field_tps = ptr->GetFieldTypes();

  int n_fields = field_tps.size();

  // convert to Rcpp::StringVector
  StringVector sv_field_tps(n_fields);

  for (int i=0; i<n_fields; ++i) {
    sv_field_tps[i] = field_tps[i];
  }

  return sv_field_tps;
}

//  [[Rcpp::export]]
NumericVector p_GeoDa__GetNumericCol(SEXP xp, std::string col_name)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<double> vals = ptr->GetNumericCol(col_name);

  int n_rows = vals.size();

  // convert to Rcpp::StringVector
  NumericVector nv_vals(n_rows);

  for (int i=0; i<n_rows; ++i) {
    nv_vals[i] = vals[i];
  }

  return nv_vals;
}

//  [[Rcpp::export]]
NumericVector p_GeoDa__GetIntegerCol(SEXP xp, std::string col_name)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<long long> vals = ptr->GetIntegerCol(col_name);

  int n_rows = vals.size();

  // convert to Rcpp::StringVector
  NumericVector nv_vals(n_rows);

  for (int i=0; i<n_rows; ++i) {
    nv_vals[i] = vals[i];
  }

  return nv_vals;
}

//  [[Rcpp::export]]
StringVector p_GeoDa__GetStringCol(SEXP xp, std::string col_name)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<std::string> vals = ptr->GetStringCol(col_name);

  int n_rows = vals.size();

  // convert to Rcpp::StringVector
  StringVector sv_vals(n_rows);

  for (int i=0; i<n_rows; ++i) {
    sv_vals[i] = vals[i];
  }

  return sv_vals;
}

// [[Rcpp::export]]
LogicalVector p_GeoDa__GetNullValues(SEXP xp, std::string col_name)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp);

  // invoke the function
  std::vector<bool> vals = ptr->GetNullValues(col_name);

  int n_rows = vals.size();

  // convert to Rcpp::LogicalVector
  LogicalVector lv_vals(n_rows);

  for (int i=0; i<n_rows; ++i) {
    lv_vals[i] = vals[i];
  }

  return lv_vals;
}

// [[Rcpp::export]]
SEXP p_GeoDa__GetPointer(SEXP xp)
{
  // return c++ object pointer
  return xp;
} 
