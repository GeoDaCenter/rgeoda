#################################################################
#' @title Weight class (Internally Used)
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
#' @field has_isolations If the weights matrix has any isolations
#' @export
Weight <- setRefClass("Weight",
  fields = list(
    gda_w = "p_GeoDaWeight",
    is_symmetric = "logical",
    sparsity = "numeric",
    density = "numeric",
    min_neighbors = "integer",
    max_neighbors = "integer",
    num_obs = "integer",
    mean_neighbors = "numeric",
    median_neighbors = "numeric",
    has_isolations = "logical"
  ),
  methods = list(
    initialize = function(o_gda_w) {
      "Constructor with a GeoDaWeight object (internally used)"
      .self$gda_w = o_gda_w
      .self$is_symmetric = gda_w$IsSymmetric()
      .self$sparsity = gda_w$GetSparsity()
      .self$density = gda_w$GetDensity()
      .self$min_neighbors = gda_w$GetMinNeighbors()
      .self$max_neighbors = gda_w$GetMaxNeighbors()
      .self$mean_neighbors = gda_w$GetMeanNeighbors()
      .self$median_neighbors = gda_w$GetMedianNeighbors()
      .self$num_obs = gda_w$GetNumObs()
      .self$has_isolations = gda_w$HasIsolations()
    },
    IsSymmetric = function() {
      "Check if weights matrix is symmetric"
      return(gda_w$IsSymmetric())
    },
    HasIsolations = function() {
      "Check if weights matrix has isolates, or if any observation has no neighbors"
      return(gda_w$HasIsolations())
    },
    GetSparsity = function() {
      "Get sparsity computed from weights matrix"
      return(gda_w$GetSparsity())
    },
    GetDensity = function() {
      "Get density computed from weights matrix"
      return (gda_w$GetDensity())
    },
    GetNeighborSize = function(idx) {
      return (gda_w$GetNeighborSize(idx))
    },
    GetNeighbors = function(idx) {
      "Get neighbors for idx-th observation, idx starts from 0"
      return (gda_w$GetNeighbors(idx))
    },
    GetNeighborWeights = function(idx) {
      "Get weights values of neighbors for idx-th observation, idx starts from 0"
      return (gda_w$GetNeighborWeights(idx))
    },
    SpatialLag = function(idx, values) {
      "Compute spatial lag values of idx-th observation, idx starts from 0"
      return (gda_w$SpatialLag(idx,values))
    },
    SaveToFile = function(out_path, layer_name, id_name, id_values) {
      "Save current spatial weights to a file.\\cr \\cr
        out_path: The path of an output weights file \\cr
        layer_name : The name of the layer of input dataset \\cr
        id_name : The id name (or field name), which is an associated column contains unique values, that makes sure that the weights are connected to the correct observations in the data table.\\cr
        id_values : The tuple of values of selected id_name (column/field)"

      return (gda_w$SaveToFile(out_path, layer_name, id_name, id_values))
    },
    GetPointer = function() {
      "Get the C++ object pointer (internally used)"
      return(gda_w$GetPointer())
    }
  )
)

############################################################
#' @title Summary of Spatial Weights
#' @description Override the summary() function for spatial weights
#' @param object A Weight object
#' @param ... summary optional parameters
#' @return A summary description of Weight object
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' summary(queen_w)
#' }
#' @export
summary.Weight <- function(object, ...) {
  gda_w <- object
  name <- c("number of observations:",
            "is symmetric: ",
            "sparsity:",
            "density:",
            "# min neighbors:",
            "# max neighbors:",
            "# mean neighbors:",
            "# median neighbors:",
            "has isolations:"
            )
  value <- c( gda_w$num_obs,
              gda_w$is_symmetric,
              gda_w$sparsity,
              gda_w$density,
              gda_w$min_neighbors,
              gda_w$max_neighbors,
              gda_w$mean_neighbors,
              gda_w$median_neighbors,
              gda_w$has_isolations)

  output <- data.frame(name, value)
  format(output)
}

#################################################################
#' @title Check if weights matrix is symmetric
#' @description Check if weights matrix is symmetric
#' @param gda_w A Weight object
#' @return Boolean value if weights matrix is symmetric
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' is_symmetric(queen_w)
#' }
#' @export
is_symmetric <- function(gda_w) {
  return(gda_w$is_symmetric)
}

#################################################################
#' @title Check if weights matrix has isolates
#' @description Check if weights matrix has isolates, or if any observation has no neighbors
#' @param gda_w A Weight object
#' @return Boolean value if weights matrix is symmetric
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' has_isolates(queen_w)
#' }
#' @export
has_isolates <- function(gda_w) {
  return(gda_w$HasIsolations())
}

#################################################################
#' @title Sparsity of a weights matrix
#' @description Get sparsity computed from weights matrix
#' @param gda_w A Weight object
#' @return Value of the weight matrix sparsity
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' weights_sparsity(queen_w)
#' }
#' @export
weights_sparsity <- function(gda_w) {
  return(gda_w$sparsity)
}

#################################################################
#' @title Density of a weights matrix
#' @description Get density computed from weights matrix
#' @param gda_w A Weight object
#' @return Value of the weight matrix density
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' weights_density(queen_w)
#' }
#' @export
weights_density <- function(gda_w) {
  return (gda_w$density)
}

#################################################################
#' @title Get neighbors for idx-th observation based on weights matrix
#' @description Get neighbors for idx-th observation, idx starts from 1
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @return Vector of the neighbor indices, idx start from 1
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' nbrs <- get_neighbors(queen_w, idx = 1)
#' cat("\nNeighbors of the 1-st observation are:", nbrs)
#' }
#' @export
get_neighbors <- function(gda_w, idx) {
  idx <- idx - 1
  nn <- gda_w$GetNeighborSize(idx)
  nbrs <- gda_w$GetNeighbors(idx)
  rtn_nbrs <- vector()
  for (i in 1:nn) {
    rtn_nbrs[i]  <- nbrs[i] + 1
  }
  return(rtn_nbrs)
}

#################################################################
#' @title Compute the spatial lag for idx-th observation and selected variable
#' @description Compute the spatial lag for idx-th observation using selected variable and current weights matrix
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @param values A vector of values
#' @return Value of the spatial lag for idx-th observation
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lag0 <- spatial_lag(queen_w, idx = 1, values = crm_prs)
#' cat("\nSpatial lag of the 1-st observation of variable crm_prs is:", lag0)
#' }
#' @export
spatial_lag <- function(gda_w, idx, values) {
  idx <- idx - 1
  return(gda_w$SpatialLag(idx, values))
}

#################################################################
#' @title Save current spatial weights to a file
#' @description Save current spatial weights to a file
#' @param gda_w A Weight object
#' @param out_path The path of an output weights file
#' @param layer_name The name of the layer of input dataset
#' @param id_name The id name (or field name), which is an associated column contains unique values, that makes sure that the weights are connected to the correct observations in the data table.
#' @param id_values The tuple of values of selected id_name (column or field)
#' @return Boolean value indicates if save successfully or failed
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' save_weights(rook_w, out_path = '/Users/xun/Downloads/Guerry_r.gal', 
#'             layer_name = 'Guerry', 
#'             id_name = 'CODE_DE', 
#'             id_values = as.integer(guerry_df['CODE_DE'][,1]))
#' }
#' @export
save_weights <- function(gda_w, out_path, layer_name, id_name, id_values) {
  return(gda_w$SaveToFile(out_path, layer_name, id_name, id_values))
}

#################################################################
#' @title Get minimum threshold of distance that makes sure each observation has at least one neighbor
#' @param geoda_obj An object of [geoda] class
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return dist A real value of minimum threshold of distance
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' dist_thres <- min_distthreshold(guerry)
#' dist_thres
#' }
#' @export
min_distthreshold <- function(geoda_obj, is_arc = FALSE, is_mile = TRUE) {
  return (p_gda_min_distthreshold(geoda_obj$GetPointer(), is_arc, is_mile))
}

#################################################################
#' @title Create a Queen contiguity weights
#' @description Create a Queen contiguity weights with options of "order", "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape file is insufficient to allow for an exact match of coordinates to determine which polygons are neighbors
#' @return w An object of Weight class
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' summary(queen_w)
#' @export
queen_weights <- function(geoda_obj, order=1, include_lower_order = FALSE, precision_threshold = 0) {

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

  gda <- geoda_obj$GetPointer()
  w <- p_gda_queen_weights(gda, order, include_lower_order, precision_threshold)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title Create a Rook contiguity weights
#' @description Create a Rook contiguity weights with options of "order", "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape file is insufficient to allow for an exact match of coordinates to determine which polygons are neighbors
#' @return w An object of Weight class
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' rook_w <- rook_weights(guerry)
#' summary(rook_w)
#' @export
rook_weights <- function(geoda_obj, order = 1, include_lower_order = FALSE, precision_threshold = 0) {

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

  w <- p_gda_rook_weights(geoda_obj$GetPointer(), order, include_lower_order, precision_threshold)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title   Create a distance-based weights
#' @description Create a distance-based weights
#' @keywords distance weights
#' @param geoda_obj An instance of geoda
#' @param dist_thres A positive numeric value of distance threshold
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value
#' @param power (optional) The power (or exponent) of a number indicates how many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return w An instance of GeoDaWeight
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' dist_thres <- min_distthreshold(guerry)
#' dist_w <- distance_weights(guerry, dist_thres)
#' summary(dist_w)
#' @export
distance_weights <- function(geoda_obj, dist_thres, power = 1.0, is_inverse = FALSE, is_arc = FALSE, is_mile=TRUE){

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (dist_thres <= 0) {
    stop("dist_thres has to be a positive numeric number.")
  }

  w <- p_gda_distance_weights(geoda_obj$GetPointer(), dist_thres, power, is_inverse, is_arc, is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
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
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' bandwidth <- min_distthreshold(guerry)
#' kernel_w <- kernel_weights(guerry, bandwidth, kernel_method = "uniform")
#' summary(kernel_w)
#' @export
kernel_weights <- function(geoda_obj, bandwidth, kernel_method,
                           use_kernel_diagonals = FALSE, power = 1.0, is_inverse = FALSE,
                           is_arc = FALSE, is_mile = TRUE) {

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

  w <- p_gda_kernel_weights(geoda_obj$GetPointer(), bandwidth, kernel_method, use_kernel_diagonals, power, is_inverse, is_arc, is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
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
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' knn6_w <- knn_weights(guerry, 6)
#' summary(knn6_w)
#' @export
knn_weights <- function(geoda_obj, k, power = 1.0, is_inverse = FALSE,
                        is_arc = FALSE, is_mile = TRUE) {

  # test if geoda_obj object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("geoda_obj object is not valid.")
  }

  if (k <= 0) {
    stop("k has to be a positive integernumeric.")
  }

  w <- p_gda_knn_weights(geoda_obj$GetPointer(), k, power, is_inverse, is_arc, is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
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
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' adptkernel_w = kernel_knn_weights(guerry, 6, "uniform")
#' summary(adptkernel_w)
#' @export
kernel_knn_weights <- function(geoda_obj, k, kernel_method, adaptive_bandwidth = TRUE,
                               use_kernel_diagonals = FALSE, power = 1.0, is_inverse = FALSE,
                               is_arc = FALSE, is_mile = TRUE ){

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

  w <- p_gda_kernel_knn_weights(geoda_obj$GetPointer(), k, power, is_inverse, is_arc, is_mile, kernel_method, 0, adaptive_bandwidth, use_kernel_diagonals)

  return(Weight$new(p_GeoDaWeight(w)))
}
