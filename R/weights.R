#################################################################
#' @title GeoDa Weight class
#' @description A wrapper class for GeoDaWeight class
#' @field gda_w An object of GeoDaWeight
#' @field is_symmetric If weights matrix is symmetric
#' @field sparsity Sparsity of weights matrix
#' @field density Density of weights matrix
#' @field min_neighbors Minimum number of neighbors
#' @field max_neighbors Maximum number of neighbors
#' @field num_obs Number of observations
#' @field mean_neighbors Mean number of neighbors
#' @field median_neighbors Median number of neighbors
#' @export
Weight <- setRefClass("Weight",
  fields = list(
    gda_w = "_p_GeoDaWeight",
    is_symmetric = "logical",
    sparsity = "numeric",
    density = "numeric",
    min_neighbors = "integer",
    max_neighbors = "integer",
    num_obs = "integer",
    mean_neighbors = "numeric",
    median_neighbors = "numeric"
  ),
  methods = list(
    initialize = function(o_gda_w) {
      "Constructor with a GeoDaWeight object (internally used)"
      .self$gda_w = o_gda_w
      .self$is_symmetric = o_gda_w$is_symmetric
      .self$sparsity = o_gda_w$sparsity
      .self$density = o_gda_w$density
      .self$min_neighbors = o_gda_w$min_nbrs
      .self$max_neighbors = o_gda_w$max_nbrs
      .self$mean_neighbors = o_gda_w$mean_nbrs
      .self$median_neighbors = o_gda_w$median_nbrs
      .self$num_obs = o_gda_w$num_obs
    },
    IsSymmetric = function() {
      "Check if weights matrix is symmetric"
      return(gda_w$is_symmetric)
    },
    HasIsolates = function() {
      "Check if weights matrix has isolates, or if any observation has no neighbors"
      return(gda_w$HasIsolates())
    },
    GetSparsity = function() {
      "Get sparsity computed from weights matrix"
      return(gda_w$sparsity)
    },
    GetDensity = function() {
      "Get density computed from weights matrix"
      return (gda_w$density)
    },
    GetNeighbors = function(idx) {
      "Get neighbors for idx-th observation, idx starts from 0"
      nn <- gda_w$GetNbrSize(idx)
      nbrs <- gda_w$GetNeighbors(idx)
      rtn_nbrs <- vector()
      for (i in 1:nn) {
        rtn_nbrs[i]  <- nbrs[i]
      }
      return(rtn_nbrs)
    },
    SpatialLag = function(idx, values) {
      "Compute spatial lag values of idx-th observation, idx starts from 0"
      return(gda_w$SpatialLag(idx,values))
    },
    SaveToFile = function(out_path, layer_name, id_name, id_values) {
      "Save current spatial weights to a file.\\cr \\cr
        out_path: The path of an output weights file \\cr
        layer_name : The name of the layer of input dataset \\cr
        id_name : The id name (or field name), which is an associated column contains unique values, that makes sure that the weights are connected to the correct observations in the data table.\\cr
        id_values : The tuple of values of selected id_name (column/field)"
      return(gda_w$Save(out_path, layer_name, id_name, id_values))
    }
  )
)

#################################################################
#' @title Get minimum threshold of distance that makes sure each observation has at least one neighbor
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return dist A real value of minimum threshold of distance
#' @export
min_distthreshold <- function(geoda_obj, ...) {
  kwargs <- list(...)
  is_arc <- ifelse(hasArg("is_arc"), kwargs$is_arc, FALSE)
  is_mile <- ifelse(hasArg("is_mile"), kwargs$is_mile, TRUE)

  return (gda_min_distthreshold(geoda_obj$gda, is_arc, is_mile))
}

#################################################################
#' @title Create a Queen contiguity weights
#' @description Create a Queen contiguity weights with options of "order", "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape file is insufficient to allow for an exact match of coordinates to determine which polygons are neighbors
#' @return w An object of Weight class
#' @export
queen_weights <- function(geoda_obj, ...) {
  kwargs <- list(...)
  order <- ifelse(hasArg("order"), kwargs$order, 1)
  include_lower_order <- ifelse(hasArg("include_lower_order"), kwargs$is_inverse, FALSE)
  precision_threshold <- ifelse(hasArg("precision_threshold"), kwargs$precision_threshold, 1.0)

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (order < 1) {
    stop("order has to be a positive integer number.")
  }

  if (precision_threshold < 0) {
    stop("precision_threshold has to be a positive numeric number.")
  }

  w <- gda_queen_weights(geoda_obj$gda, "", order, include_lower_order, precision_threshold)
  return(Weight$new(w))
}

#################################################################
#' @title Create a Rook contiguity weights
#' @description Create a Rook contiguity weights with options of "order", "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape file is insufficient to allow for an exact match of coordinates to determine which polygons are neighbors
#' @return w An object of Weight class
#' @export
rook_weights <- function(geoda_obj, ...) {
  kwargs <- list(...)
  order <- ifelse(hasArg("order"), kwargs$order, 1)
  include_lower_order <- ifelse(hasArg("include_lower_order"), kwargs$is_inverse, FALSE)
  precision_threshold <- ifelse(hasArg("precision_threshold"), kwargs$precision_threshold, 1.0)

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (order < 1) {
    stop("order has to be a positive integer number.")
  }

  if (precision_threshold < 0) {
    stop("precision_threshold has to be a positive numeric number.")
  }

  w <- gda_rook_weights(geoda_obj$gda, "", order, include_lower_order, precision_threshold)
  return(Weight$new(w))
}

#################################################################
#' @title   Create a distance-based weights
#' @description Create a distance-based weights
#' @keywords distance weights
#' @param geoda_obj An instance of geoda
#' @param dist_thres A positive numeric value of distance threshold
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value
#' @param power (optional) The power (or exponent) of a number says how many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return w An instance of GeoDaWeight
#' @export
distance_weights <- function(geoda_obj, dist_thres, ...) {
  kwargs <- list(...)
  power <- ifelse(hasArg("power"), kwargs$power, 1.0)
  is_inverse <- ifelse(hasArg("is_inverse"), kwargs$is_inverse, FALSE)
  is_arc <- ifelse(hasArg("is_arc"), kwargs$is_arc, FALSE)
  is_mile <- ifelse(hasArg("is_mile"), kwargs$is_mile, TRUE)

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (dist_thres <= 0) {
    stop("dist_thres has to be a positive numeric number.")
  }

  w <- gda_distance_weights(geoda_obj$gda, dist_thres, "", power, is_inverse, is_arc, is_mile, "", FALSE)
  return(Weight$new(w))
}

#################################################################
#' @title   Create a kernel weights with fixed bandwidth
#' @description Create a kernel weights by specifying a bandwidth and a kernel method
#' @keywords kernel weights
#' @param geoda_obj An instance of GeoDa
#' @param bandwidth A positive numeric value of bandwidth
#' @param kernel_method a string value, which has to be one of 'triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value
#' @param power (optional) The power (or exponent) of a number says how many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return w An instance of Weight
#' @export
kernel_weights <- function(geoda_obj, bandwidth, kernel_method, ...) {
  kwargs <- list(...)
  use_kernel_diagonals <- ifelse(hasArg("use_kernel_diagonals"), kwargs$use_kernel_diagonals, FALSE)
  power <- ifelse(hasArg("power"), kwargs$power, 1.0)
  is_inverse <- ifelse(hasArg("is_inverse"), kwargs$is_inverse, FALSE)
  is_arc <- ifelse(hasArg("is_arc"), kwargs$is_arc, FALSE)
  is_mile <- ifelse(hasArg("is_mile"), kwargs$is_mile, TRUE)

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (bandwidth <= 0) {
    stop("bandwidth has to be a positive numeric number.")
  }

  if (!(kernel_method %in% c('triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'))) {
    stop("kernel_method has to be one of 'triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'.")
  }

  w <- gda_distance_weights(geoda_obj$gda, bandwidth, "", power, is_inverse, is_arc, is_mile, kernel_method, use_kernel_diagonals)
  return(Weight$new(w))
}

#################################################################
#' @title  Create a KNN weights
#' @description Create a k-nearest neighbors based weights
#' @keywords weights knn
#' @param geoda_obj An instance of geoda
#' @param k a positive integer number for k-nearest neighbors
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value
#' @param power (optional) The power (or exponent) of a number says how many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return w An instance of Weight
#' @export
knn_weights <- function(geoda_obj, k, ...) {
  kwargs <- list(...)
  power <- ifelse(hasArg("power"), kwargs$power, 1.0)
  is_inverse <- ifelse(hasArg("is_inverse"), kwargs$is_inverse, FALSE)
  is_arc <- ifelse(hasArg("is_arc"), kwargs$is_arc, FALSE)
  is_mile <- ifelse(hasArg("is_mile"), kwargs$is_mile, TRUE)

  # test if geoda_obj object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("geoda_obj object is not valid.")
  }

  if (k <= 0) {
    stop("k has to be a positive integernumeric.")
  }

  w <- gda_knn_weights(geoda_obj$gda, k, power, is_inverse, is_arc, is_mile, "", 0, FALSE, FALSE, "")
  return(Weight$new(w))
}

#################################################################
#' @title  Create a kernel weights with KNN specified bandwidth
#' @description Create a kernel weights by specifying k-nearest neighbors and a kernel method
#' @keywords kernel weights knn
#' @param geoda_obj An instance of geoda
#' @param k a positive integer number for k-nearest neighbors
#' @param kernel_method a string value, which has to be one of 'triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param adaptive_bandwidth (optional) TRUE (default) or FALSE: TRUE use adaptive bandwidth calculated using distance of k-nearest neithbors,
#' FALSE use max distance of all observation to their k-nearest neighbors
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value
#' @param power (optional) The power (or exponent) of a number says how many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return w An instance of Weight
#' @export
kernel_knn_weights <- function(geoda_obj, k, kernel_method, ...) {
  kwargs <- list(...)
  adaptive_bandwidth <- ifelse(hasArg("adaptive_bandwidth"), kwargs$adaptive_bandwidth, TRUE)
  use_kernel_diagonals <- ifelse(hasArg("use_kernel_diagonals"), kwargs$use_kernel_diagonals, FALSE)
  power <- ifelse(hasArg("power"), kwargs$power, 1.0)
  is_inverse <- ifelse(hasArg("is_inverse"), kwargs$is_inverse, FALSE)
  is_arc <- ifelse(hasArg("is_arc"), kwargs$is_arc, FALSE)
  is_mile <- ifelse(hasArg("is_mile"), kwargs$is_mile, TRUE)

  # test if geoda_obj object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("geoda_obj object is not valid.")
  }

  if (k <= 0) {
    stop("k has to be a positive integernumeric.")
  }

  if (!(kernel_method %in% c('triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'))) {
    stop("kernel_method has to be one of 'triangular', 'uniform', 'epanechnikov', 'quartic', 'gaussian'.")
  }

  w <- gda_knn_weights(geoda_obj$gda, k, power, is_inverse, is_arc, is_mile, kernel_method, 0, adaptive_bandwidth, use_kernel_diagonals, "")
  return(Weight$new(w))
}
