// This file is used to wrap C++ classes and functions defines in RcppExports.R
// All other R script files will use this file as a bridge to C++ classes and functions
//
// Author: lixun910@gmail.com
// Changes:
// 10/29/2020 init rcpp_weights.cpp

#include <Rcpp.h>
#include "libgeoda_src/weights/GeodaWeight.h"
#include "libgeoda_src/gda_weights.h"
#include "libgeoda_src/libgeoda.h"

using namespace Rcpp;

// [[Rcpp::export]]
SEXP p_GeoDaWeight__GetPointer(SEXP xp)
{
  // return c++ object pointer
  return xp;
} 

//  [[Rcpp::export]]
int p_GeoDaWeight__GetNumObs(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  return ptr->GetNumObs();
}

//  [[Rcpp::export]]
bool p_GeoDaWeight__IsSymmetric(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  bool is_sym = ptr->IsSymmetric();

  return is_sym;
}

//  [[Rcpp::export]]
bool p_GeoDaWeight__HasIsolations(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  bool has_iso = ptr->HasIsolations();

  return has_iso;
}

//  [[Rcpp::export]]
double p_GeoDaWeight__GetSparsity(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  double sparsity = ptr->GetSparsity();

  return sparsity;
}

//  [[Rcpp::export]]
double p_GeoDaWeight__GetDensity(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  double density = ptr->GetDensity();

  return density;
}

//  [[Rcpp::export]]
int p_GeoDaWeight__GetMinNeighbors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  int min_nbrs = ptr->GetMinNbrs();

  return min_nbrs;
}

//  [[Rcpp::export]]
int p_GeoDaWeight__GetMaxNeighbors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  int max_nbrs = ptr->GetMaxNbrs();

  return max_nbrs;
}

//  [[Rcpp::export]]
double p_GeoDaWeight__GetMeanNeighbors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  return ptr->GetMeanNbrs();
}

//  [[Rcpp::export]]
double p_GeoDaWeight__GetMedianNeighbors(SEXP xp)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  return ptr->GetMedianNbrs();
}

//  [[Rcpp::export]]
double p_GeoDaWeight__SpatialLag(SEXP xp, int obs_idx, SEXP vals)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  std::vector<double> _vals = as<std::vector<double> >(vals);

  // invoke the function
  double lag = ptr->SpatialLag(obs_idx, _vals);

  return lag;
}

//  [[Rcpp::export]]
int p_GeoDaWeight__GetNeighborSize(SEXP xp, int obs_idx)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  int nn_sz = ptr->GetNbrSize(obs_idx);

  return nn_sz;
}

//  [[Rcpp::export]]
bool p_GeoDaWeight__SaveToFile(SEXP xp, std::string out_path, std::string layer_name, std::string id_name, SEXP id_vec)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  if (TYPEOF(id_vec) == INTSXP) { // using integers as id_vec
    // convert
    std::vector<int> _id_vec = as<std::vector<int> >(id_vec);

    // invoke the function
    return ptr->Save(out_path.c_str(), layer_name.c_str(), id_name.c_str(), _id_vec);
  } else {
    // using strings as id_vec
    std::vector<std::string> tmp = as<std::vector<std::string> >(id_vec);
    std::vector<const char*> _id_vec;
    for( int i=0; i < tmp.size(); i++ ){
      _id_vec.push_back(tmp[i].c_str());
    }

    // invoke the function
    return ptr->Save(out_path.c_str(), layer_name.c_str(), id_name.c_str(), _id_vec);
  }
}

//  [[Rcpp::export]]
NumericVector p_GeoDaWeight__GetNeighbors(SEXP xp, int obs_idx)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  std::vector<long> nn = ptr->GetNeighbors(obs_idx);

  // convert to Rcpp::StringVector
  NumericVector nv_nn(nn.size());

  for (int i=0; i<nn.size(); ++i) {
    nv_nn[i] = nn[i];
  }

  return nv_nn;
}

//  [[Rcpp::export]]
NumericVector p_GeoDaWeight__GetNeighborWeights(SEXP xp, int obs_idx)
{
  // grab the object as a XPtr (smart pointer) to GeoDaWeight
  Rcpp::XPtr<GeoDaWeight> ptr(xp);

  // invoke the function
  std::vector<double> nnw = ptr->GetNeighborWeights(obs_idx);

  // convert to Rcpp::StringVector
  NumericVector nv_nnw(nnw.size());

  for (int i=0; i<nnw.size(); ++i) {
    nv_nnw[i] = nnw[i];
  }

  return nv_nnw;
}

//  [[Rcpp::export]]
double p_gda_min_distthreshold(SEXP xp_geoda, bool is_arc, bool is_mile)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  double min_dist = gda_min_distthreshold(geoda, is_arc, is_mile);

  return min_dist;
}

//  [[Rcpp::export]]
SEXP p_gda_queen_weights(SEXP xp_geoda, int order, bool include_lower_order, double precision_threshold)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_queen_weights(geoda, order, include_lower_order, precision_threshold);

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true); // true: we need to register a delete finalizer with the external pointer.

  return w_ptr;
}

//  [[Rcpp::export]]
SEXP p_gda_rook_weights(SEXP xp_geoda, int order, bool include_lower_order, double precision_threshold)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_rook_weights(geoda, order, include_lower_order, precision_threshold);

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true);
  return w_ptr;
}

//  [[Rcpp::export]]
SEXP p_gda_distance_weights(SEXP xp_geoda, double dist_thres, double power, bool is_inverse, bool is_arc, bool is_mile)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_distance_weights(geoda, dist_thres, "", power, is_inverse, is_arc, is_mile, "", FALSE);

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true);
  return w_ptr;
}

//  [[Rcpp::export]]
SEXP p_gda_kernel_weights(SEXP xp_geoda, double bandwidth, std::string kernel_method, bool use_kernel_diagonals, double power, bool is_inverse, bool is_arc, bool is_mile)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_distance_weights(geoda, bandwidth, "", power, is_inverse, is_arc, is_mile, kernel_method, use_kernel_diagonals);

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true);
  return w_ptr;
}

//  [[Rcpp::export]]
SEXP p_gda_knn_weights(SEXP xp_geoda, int k, double power, bool is_inverse, bool is_arc, bool is_mile)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_knn_weights(geoda, k, power, is_inverse, is_arc, is_mile, "", 0, FALSE, FALSE, "");

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true);
  return w_ptr;
}

//  [[Rcpp::export]]
SEXP p_gda_kernel_knn_weights(SEXP xp_geoda, int k, double power, bool is_inverse, bool is_arc, bool is_mile, std::string kernel_method, double bandwidth, bool adaptive_bandwidth, bool use_kernel_diagonals)
{
  // grab the object as a XPtr (smart pointer) to GeoDa
  Rcpp::XPtr<GeoDa> ptr(xp_geoda);
  GeoDa* geoda = static_cast<GeoDa*> (R_ExternalPtrAddr(ptr));

  // invoke the function
  GeoDaWeight* w = gda_knn_weights(geoda, k, power, is_inverse, is_arc, is_mile, kernel_method, 0, adaptive_bandwidth, use_kernel_diagonals, "");

  Rcpp::XPtr<GeoDaWeight> w_ptr(w, true);
  return w_ptr;
}
