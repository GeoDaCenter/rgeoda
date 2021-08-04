#################################################################
#' @title Weight class (Internally Used)
#' @description A wrapper class for p_GeoDaWeight class
#' @field gda_w An object of p_GeoDaWeight-class
#' @field is_symmetric If weights matrix is symmetric
#' @field sparsity Sparsity of weights matrix
#' @field min_neighbors Minimum number of neighbors
#' @field max_neighbors Maximum number of neighbors
#' @field num_obs Number of observations
#' @field mean_neighbors Mean number of neighbors
#' @field median_neighbors Median number of neighbors
#' @field has_isolates If the weights matrix has any isolates
#' @export
Weight <- setRefClass("Weight",
  fields = list(
    gda_w = "p_GeoDaWeight",
    is_symmetric = "logical",
    sparsity = "numeric",
    min_neighbors = "integer",
    max_neighbors = "integer",
    num_obs = "integer",
    mean_neighbors = "numeric",
    median_neighbors = "numeric",
    has_isolates= "logical"
  ),
  methods = list(
    initialize = function(o_gda_w) {
      "Constructor with a GeoDaWeight object (internally used)"
      .self$gda_w = o_gda_w
      .self$Update(FALSE)
    },
    Update = function(updateStats = TRUE) {
      "Update the weights meta data"
      if (updateStats == TRUE) {
        gda_w$GetNbrStats()
      }
      .self$is_symmetric = gda_w$IsSymmetric()
      .self$sparsity = gda_w$GetSparsity()
      .self$min_neighbors = gda_w$GetMinNeighbors()
      .self$max_neighbors = gda_w$GetMaxNeighbors()
      .self$mean_neighbors = gda_w$GetMeanNeighbors()
      .self$median_neighbors = gda_w$GetMedianNeighbors()
      .self$num_obs = gda_w$GetNumObs()
      .self$has_isolates = gda_w$HasIsolates()
    },
    SetNeighbors = function(idx, nbrs) {
      "Set neighbors for one observation"
      gda_w$SetNeighbors(idx, nbrs)
    },
    SetNeighborsAndWeights = function(idx, nbrs, nbr_w) {
      "Set neighbors with weights values for one observation"
      gda_w$SetNeighborsAndWeights(idx, nbrs, nbr_w)
    },
    IsSymmetric = function() {
      "Check if weights matrix is symmetric"
      return(gda_w$IsSymmetric())
    },
    HasIsolates = function() {
      "Check if weights matrix has isolates, or if any observation has no
      neighbors"
      return(gda_w$HasIsolates())
    },
    GetSparsity = function() {
      "Get sparsity computed from weights matrix"
      return(gda_w$GetSparsity())
    },
    GetNeighborSize = function(idx) {
      return (gda_w$GetNeighborSize(idx))
    },
    GetNeighbors = function(idx) {
      "Get neighbors for idx-th observation, idx starts from 0"
      return (gda_w$GetNeighbors(idx))
    },
    GetNeighborWeights = function(idx) {
      "Get weights values of neighbors for idx-th observation,
      idx starts from 0"
      return (gda_w$GetNeighborWeights(idx))
    },
    SpatialLag = function(values) {
      "Compute spatial lag values for values of selected variable"
      return (gda_w$SpatialLag(values))
    },
    SaveToFile = function(out_path, layer_name, id_name, id_values) {
      "Save current spatial weights to a file.\\cr \\cr
        out_path: The path of an output weights file \\cr
        layer_name : The name of the layer of input dataset \\cr
        id_name : The id name (or field name), which is an associated column
        contains unique values, that makes sure that the weights are connected
        to the correct observations in the data table.\\cr
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
#' @return A summary description of an instance of Weight-class
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' summary(queen_w)
#' }
#' @export
summary.Weight <- function(object, ...) {
  gda_w <- object
  name <- c("number of observations:",
            "is symmetric: ",
            "sparsity:",
            "# min neighbors:",
            "# max neighbors:",
            "# mean neighbors:",
            "# median neighbors:",
            "has isolates:"
            )
  value <- c( as.character(gda_w$num_obs),
              as.character(gda_w$is_symmetric),
              as.character(gda_w$sparsity),
              as.character(gda_w$min_neighbors),
              as.character(gda_w$max_neighbors),
              as.character(gda_w$mean_neighbors),
              as.character(gda_w$median_neighbors),
              as.character(gda_w$has_isolates))

  output <- data.frame(name, value)
  format(output)
}

#################################################################
#' @title Create an empty weights
#' @description Create an empty weights
#' @param num_obs The number of observations for this empty weights
#' @return An instance of Weight-class
#' @export
create_weights <- function(num_obs) {
  return (Weight$new(p_GeoDaWeight(num_obs)))
}

#################################################################
#' @title Set neighbors of an observation
#' @description Set neighbors for idx-th observation, idx starts from 1
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @param nbrs A list indicates the neighbors of idx-th observation
#' (id start from 1)
#' @examples
#' \dontrun{
#' new_w <- create_weights(10)
#' set_neighbors(new_w, 1, c(2,3))
#' update_weights(new_w)
#' }
#' @export
set_neighbors <- function(gda_w, idx, nbrs) {
  gda_w$SetNeighbors(idx, nbrs)
}

#################################################################
#' @title Set neighbors and weights values of an observation
#' @description Set neighbors and the associated weights values for idx-th
#' observation, idx starts from 1
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @param nbrs A list indicates the neighbors of idx-th observation
#' (id start from 1)
#' @param wvals A list indicates the associated weights values of the neighbors
#' @examples
#' \dontrun{
#' new_w <- create_weights(10)
#' set_neighbors(new_w, 1, c(2,3))
#' update_weights(new_w)
#' }
#' @export
set_neighbors_with_weights <- function(gda_w, idx, nbrs, wvals) {
  gda_w$SetNeighborsAndWeights(idx, nbrs, wvals)
}

#################################################################
#' @title Update meta data of a spatial weights
#' @description Update meta data of a spatial weights. This function can be used
#' after calling `set_neighbor()` function .
#' @param gda_w A Weight object
#' @examples
#' \dontrun{
#' new_w <- create_weights(10)
#' set_neighbors(new_w, 1, c(2,3))
#' update_weights(new_w)
#' }
#' @export
update_weights <- function(gda_w) {
  gda_w$Update()
}

#################################################################
#' @title Symmetry of Weights Matrix
#' @description Check if weights matrix is symmetric
#' @param gda_w A Weight object
#' @return A boolean value indicates if weights matrix is symmetric
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' is_symmetric(queen_w)
#' }
#' @export
is_symmetric <- function(gda_w) {
  return(gda_w$is_symmetric)
}

#################################################################
#' @title Isolation/Island in Spatial Weights
#' @description Check if weights matrix has isolates, or if any observation has
#' no neighbors
#' @param gda_w A Weight object
#' @return A boolean value indicates if weights matrix is symmetric
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' has_isolates(queen_w)
#' }
#' @export
has_isolates <- function(gda_w) {
  return(gda_w$HasIsolates())
}

#################################################################
#' @title Sparsity of Spatial Weights
#' @description Get sparsity (% Non-zero) computed from weights matrix
#' @param gda_w A Weight object
#' @return A numeric value of spatial weights sparsity
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' weights_sparsity(queen_w)
#' }
#' @export
weights_sparsity <- function(gda_w) {
  return(gda_w$sparsity)
}

#################################################################
#' @title Neighbors of one observation
#' @description Get neighbors for idx-th observation, idx starts from 1
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @return A numeric vector of the neighbor indices, which start from 1
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
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
#' @title Weights values of the neighbors of one observation
#' @description Get the associated weights values of neighbors for idx-th
#' observation
#' @param gda_w A Weight object
#' @param idx A value indicates idx-th observation, idx start from 1
#' @return A numeric vector of the neighbor indices, which start from 1
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' nbrs <- get_neighbors_weights(queen_w, idx = 1)
#' cat("\nNeighbors of the 1-st observation are:", nbrs)
#' }
#' @export
get_neighbors_weights <- function(gda_w, idx) {
  idx <- idx - 1
  nn <- gda_w$GetNeighborSize(idx)
  nbrs <- gda_w$GetNeighborWeights(idx)
  rtn_nbrs <- vector()
  for (i in 1:nn) {
    rtn_nbrs[i]  <- nbrs[i]
  }
  return(rtn_nbrs)
}

#################################################################
#' @title Spatial Lag
#' @description Compute the spatial lag for idx-th observation using selected
#' variable and current weights matrix
#' @param gda_w A Weight object
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @return A data.frame with one column "Spatial Lag"
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' crm_lag <- spatial_lag(queen_w, guerry["Crm_prs"])
#' crm_lag
#' }
#' @export
spatial_lag <- function(gda_w, df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]
  return(gda_w$SpatialLag(data))
}

#################################################################
#' @title Maximum Neighbors of Spatial Weights
#' @description Get the number of maximum neighbors of spatial weights
#' @param gda_w A Weight object
#' @return The number of maximum neighbors of spatial weights
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' max_neighbors(queen_w)
#' }
#' @export
max_neighbors <- function(gda_w) {
  return(gda_w$max_neighbors)
}

#################################################################
#' @title Minimum Neighbors of Spatial Weights
#' @description Get the number of minimum neighbors of spatial weights
#' @param gda_w A Weight object
#' @return The number of minimum neighbors of spatial weights
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' min_neighbors(queen_w)
#' }
#' @export
min_neighbors <- function(gda_w) {
  return(gda_w$min_neighbors)
}

#################################################################
#' @title Mean Neighbors of Spatial Weights
#' @description Get the number of mean neighbors of spatial weights
#' @param gda_w A Weight object
#' @return The number of mean neighbors of spatial weights
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' mean_neighbors(queen_w)
#' }
#' @export
mean_neighbors <- function(gda_w) {
  return(gda_w$mean_neighbors)
}

#################################################################
#' @title Median Neighbors of Spatial Weights
#' @description Get the number of median neighbors of spatial weights
#' @param gda_w A Weight object
#' @return The number of median neighbors of spatial weights
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' median_neighbors(queen_w)
#' }
#' @export
median_neighbors <- function(gda_w) {
  return(gda_w$median_neighbors)
}

#################################################################
#' @title Save Spatial Weights
#' @description Save spatial weights to a file
#' @param gda_w A Weight object
#' @param out_path The path of an output weights file
#' @param id_variable The id variable (a data.frame) that defines the unique
#' value of each observation when saving a weights file
#' @param layer_name (optional) The name of the layer of input dataset
#' @return A boolean value indicates if save successfully or failed
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' save_weights(quen_w, guerry_df['CODE_DE'], out_path = '/path/Guerry_r.gal')
#' }
#' @export
save_weights <- function(gda_w, id_variable, out_path, layer_name="") {
  if (inherits(id_variable, "data.frame") == FALSE) {
    stop("The id_variable needs to be a data.frame.")
  }
  id_values <- id_variable[[1]]
  id_name <- names(id_variable)[[1]]

  if (id_name == "") {
    stop("The id_variable doesn't have a column name.")
  }

  return(gda_w$SaveToFile(out_path, layer_name, id_name, id_values))
}

#################################################################
#' @title Queen Contiguity Spatial Weights
#' @description Create a Queen contiguity weights with options of "order",
#' "include lower order" and "precision threshold"
#' @param sf_obj An sf (simple feature) object
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order
#' neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape
#' file is insufficient to allow for an exact match of coordinates to determine
#' which polygons are neighbors
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' summary(queen_w)
#' @export
queen_weights <- function(sf_obj, order=1, include_lower_order = FALSE,
                          precision_threshold = 0) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_queen_weights(geoda_obj, order, include_lower_order,
                            precision_threshold))
}


#################################################################
#' @title Rook Contiguity Spatial Weights
#' @description Create a Rook contiguity weights with options of "order",
#' "include lower order" and "precision threshold"
#' @param sf_obj An sf (simple feature) object
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order
#' neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape
#' file is insufficient to allow for an exact match of coordinates to determine
#' which polygons are neighbors
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' rook_w <- rook_weights(guerry)
#' summary(rook_w)
#' @export
rook_weights <- function(sf_obj, order = 1, include_lower_order = FALSE,
                         precision_threshold = 0) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_rook_weights(geoda_obj, order, include_lower_order,
                           precision_threshold))
}

#################################################################
#' @title Minimum Distance Threshold for Distance-based Weights
#' @description Get minimum threshold of distance that makes sure each
#' observation has at least one neighbor
#' @param sf_obj An sf (simple feature) object
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, if 'is_arc' option is
#' TRUE, then 'is_mile' will set distance unit to 'mile' or 'km'.
#' @return A numeric value of minimum threshold of distance
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' dist_thres <- min_distthreshold(guerry)
#' dist_thres
#' }
#' @export
min_distthreshold <- function(sf_obj, is_arc = FALSE, is_mile = TRUE) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_min_distthreshold(geoda_obj, is_arc, is_mile))
}

#################################################################
#' @title Distance-based Spatial Weights
#' @description Create a distance-based weights
#' @keywords distance weights
#' @param sf_obj An sf (simple feature) object
#' @param dist_thres A positive numeric value of distance threshold
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' dist_thres <- min_distthreshold(guerry)
#' dist_w <- distance_weights(guerry, dist_thres)
#' summary(dist_w)
#' @export
distance_weights <- function(sf_obj, dist_thres, power = 1.0,
                             is_inverse = FALSE, is_arc = FALSE, is_mile=TRUE) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_distance_weights(geoda_obj, dist_thres, power, is_inverse, is_arc,
                               is_mile))
}

#################################################################
#' @title Distance-based Kernel Spatial Weights
#' @description Create a kernel weights by specifying a bandwidth and a kernel
#' method
#' @keywords kernel weights
#' @param sf_obj An sf (simple feature) object
#' @param bandwidth A positive numeric value of bandwidth
#' @param kernel_method a string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#' on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' bandwidth <- min_distthreshold(guerry)
#' kernel_w <- kernel_weights(guerry, bandwidth, kernel_method = "uniform")
#' summary(kernel_w)
#' @export
kernel_weights <- function(sf_obj, bandwidth, kernel_method,
                           use_kernel_diagonals = FALSE, power = 1.0,
                           is_inverse = FALSE, is_arc = FALSE, is_mile = TRUE) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_kernel_weights(geoda_obj, bandwidth, kernel_method,
                             use_kernel_diagonals, power, is_inverse, is_arc,
                             is_mile))
}

#################################################################
#' @title K-Nearest Neighbors-based Spatial Weights
#' @description Create a k-nearest neighbors based spatial weights
#' @keywords weights knn
#' @param sf_obj An sf (simple feature) object
#' @param k a positive integer number for k-nearest neighbors
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' knn6_w <- knn_weights(guerry, 6)
#' summary(knn6_w)
#' @export
knn_weights <- function(sf_obj, k, power = 1.0, is_inverse = FALSE,
                        is_arc = FALSE, is_mile = TRUE) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_knn_weights(geoda_obj, k, power, is_inverse, is_arc, is_mile))
}

#################################################################
#' @title K-NN Kernel Spatial Weights
#' @description Create a kernel weights by specifying k-nearest neighbors and a
#' kernel method
#' @keywords kernel weights knn
#' @param sf_obj An sf (simple feature) object
#' @param k a positive integer number for k-nearest neighbors
#' @param kernel_method a string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param adaptive_bandwidth (optional) TRUE (default) or FALSE: TRUE use
#' adaptive bandwidth calculated using distance of k-nearest neithbors,
#' FALSE use max distance of all observation to their k-nearest neighbors
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#'  on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' adptkernel_w = kernel_knn_weights(guerry, 6, "uniform")
#' summary(adptkernel_w)
#' @export
kernel_knn_weights <- function(sf_obj, k, kernel_method,
                               adaptive_bandwidth = TRUE,
                               use_kernel_diagonals = FALSE, power = 1.0,
                               is_inverse = FALSE,
                               is_arc = FALSE, is_mile = TRUE ) {
  geoda_obj <- getGeoDaObj(sf_obj)
  return (gda_kernel_knn_weights(geoda_obj, k, kernel_method,
                                 adaptive_bandwidth, use_kernel_diagonals,
                                 power, is_inverse, is_arc, is_mile))
}

#################################################################
# THE FOLLOWING FUNCTIONS ARE FOR INTERNALLY USE AND TEST ONLY
#################################################################
#' @title (For internally use and test only) Queen Contiguity Spatial Weights
#' @description Create a Queen contiguity weights with options of "order",
#' "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order
#' neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape
#' file is insufficient to allow for an exact match of coordinates to determine
#' which polygons are neighbors
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- gda_queen_weights(guerry)
#' summary(queen_w)
#' }
#' @export
gda_queen_weights <- function(geoda_obj, order=1, include_lower_order = FALSE,
                              precision_threshold = 0) {

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
#' @title (For internally use and test only) Rook Contiguity Spatial Weights
#' @description Create a Rook contiguity weights with options of "order",
#' "include lower order" and "precision threshold"
#' @param geoda_obj An object of [geoda] class
#' @param order  (Optional) Order of contiguity
#' @param include_lower_order (Optional)  Whether or not the lower order
#' neighbors should be included in the weights structure
#' @param precision_threshold  (Optional) The precision of the underlying shape
#'  file is insufficient to allow for an exact match of coordinates to determine
#'   which polygons are neighbors
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' rook_w <- gda_rook_weights(guerry)
#' summary(rook_w)
#' }
#' @export
gda_rook_weights <- function(geoda_obj, order = 1, include_lower_order = FALSE,
                             precision_threshold = 0) {

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

  w <- p_gda_rook_weights(geoda_obj$GetPointer(), order, include_lower_order,
                          precision_threshold)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title (For internally use and test only) Minimum Distance Threshold for
#' Distance-based Weights
#' @description Get minimum threshold of distance that makes sure each
#' observation has at least one neighbor
#' @param geoda_obj An instance of geoda-class
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, if 'is_arc' option is
#' TRUE, then 'is_mile' will set distance unit to 'mile' or 'km'.
#' @return A numeric value of minimum threshold of distance
#' @export
gda_min_distthreshold <- function(geoda_obj, is_arc = FALSE, is_mile = TRUE) {
  return (p_gda_min_distthreshold(geoda_obj$GetPointer(), is_arc, is_mile))
}


#################################################################
#' @title (For internally use and test only) Distance-based Spatial Weights
#' @description Create a distance-based weights
#' @keywords distance weights
#' @param geoda_obj An instance of geoda-class
#' @param dist_thres A positive numeric value of distance threshold
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' dist_thres <- gda_min_distthreshold(guerry)
#' dist_w <- gda_distance_weights(guerry, dist_thres)
#' summary(dist_w)
#' }
#' @export
gda_distance_weights <- function(geoda_obj, dist_thres, power = 1.0,
                                 is_inverse = FALSE, is_arc = FALSE,
                                 is_mile=TRUE) {

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (dist_thres <= 0) {
    stop("dist_thres has to be a positive numeric number.")
  }

  w <- p_gda_distance_weights(geoda_obj$GetPointer(), dist_thres, power,
                              is_inverse, is_arc, is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title (For internally use and test only) Distance-based Kernel Spatial
#' Weights
#' @description Create a kernel weights by specifying a bandwidth and a kernel
#' method
#' @keywords kernel weights
#' @param geoda_obj An instance of geoda-class
#' @param bandwidth A positive numeric value of bandwidth
#' @param kernel_method a string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#' on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' bandwidth <- gda_min_distthreshold(guerry)
#' kernel_w <- gda_kernel_weights(guerry, bandwidth, kernel_method = "uniform")
#' summary(kernel_w)
#' }
#' @export
gda_kernel_weights <- function(geoda_obj, bandwidth, kernel_method,
                           use_kernel_diagonals = FALSE, power = 1.0,
                           is_inverse = FALSE, is_arc = FALSE, is_mile = TRUE) {

  # test if gda object works
  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  if (bandwidth <= 0) {
    stop("bandwidth has to be a positive numeric number.")
  }

  if (!(kernel_method %in% c('triangular', 'uniform', 'epanechnikov',
                             'quartic', 'gaussian'))) {
    stop("kernel_method has to be one of 'triangular', 'uniform',
         'epanechnikov', 'quartic', 'gaussian'.")
  }

  w <- p_gda_kernel_weights(geoda_obj$GetPointer(), bandwidth, kernel_method,
                            use_kernel_diagonals, power, is_inverse, is_arc,
                            is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title (For internally use and test only) K-Nearest Neighbors-based Spatial
#' Weights
#' @description Create a k-nearest neighbors based spatial weights
#' @keywords weights knn
#' @param geoda_obj An instance of geoda
#' @param k a positive integer number for k-nearest neighbors
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit
#' from mile to km.
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' knn6_w <- gda_knn_weights(guerry, 6)
#' summary(knn6_w)
#' }
#' @export
gda_knn_weights <- function(geoda_obj, k, power = 1.0, is_inverse = FALSE,
                        is_arc = FALSE, is_mile = TRUE) {

  # test if geoda_obj object works
  if (geoda_obj$GetNumObs() <= 0) {
    stop("geoda_obj object is not valid.")
  }

  if (k <= 0) {
    stop("k has to be a positive integernumeric.")
  }

  w <- p_gda_knn_weights(geoda_obj$GetPointer(), k, power, is_inverse, is_arc,
                         is_mile)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title (For internally use and test only) K-NN Kernel Spatial Weights
#' @description Create a kernel weights by specifying k-nearest neighbors and a
#' kernel method
#' @keywords kernel weights knn
#' @param geoda_obj An instance of geoda
#' @param k a positive integer number for k-nearest neighbors
#' @param kernel_method a string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'
#' @param adaptive_bandwidth (optional) TRUE (default) or FALSE: TRUE use
#' adaptive bandwidth calculated using distance of k-nearest neithbors,
#' FALSE use max distance of all observation to their k-nearest neighbors
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#' on the diagonal of weights matrix
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value
#' @param power (optional) The power (or exponent) of a number says how many
#' times to use the number in a multiplication.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance
#' between two observations
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from
#'  mile to km.
#' @return An instance of Weight-class
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' adptkernel_w = gda_kernel_knn_weights(guerry, 6, "uniform")
#' summary(adptkernel_w)
#' }
#' @export
gda_kernel_knn_weights <- function(geoda_obj, k, kernel_method,
                                   adaptive_bandwidth = TRUE,
                                   use_kernel_diagonals = FALSE,
                                   power = 1.0, is_inverse = FALSE,
                                   is_arc = FALSE, is_mile = TRUE) {

  # test if geoda_obj object works
  if (geoda_obj$GetNumObs() <= 0) {
    stop("geoda_obj object is not valid.")
  }

  if (k <= 0) {
    stop("k has to be a positive integernumeric.")
  }

  if (!(kernel_method %in% c("triangular", "uniform", "epanechnikov", "quartic",
                             "gaussian"))) {
    stop("kernel_method has to be one of 'triangular', 'uniform',
         'epanechnikov', 'quartic', 'gaussian'.")
  }

  w <- p_gda_kernel_knn_weights(geoda_obj$GetPointer(), k, power, is_inverse,
                                is_arc, is_mile, kernel_method, 0,
                                adaptive_bandwidth, use_kernel_diagonals)

  return(Weight$new(p_GeoDaWeight(w)))
}


#################################################################
#' @title spatial weights to matrix
#' @description Convert a GeoDa spatial weights object to a Matrix object
#' @param x A weights object
#' @param rownames optional, a single column name or column number to use as the
#'  rownames in the returned matrix. If TRUE the key of the data.table will be
#'  used if it is a single column, otherwise the first column in the data.table
#'  will be used.
#' @param rownames.value optional, a vector of values to be used as the rownames
#'  in the returned matrix. It must be the same length as nrow(x).
#' @param ... Required to be present because the generic `as.matrix` generic
#' has it. Arguments here are not currently used or passed on by this method.
#' @return A matrix object
#' @export
as.matrix.Weight <- function(x, rownames=NULL, rownames.value=NULL, ...) {
  if (length(class(x)) == 1 && class(x) == "Weight") {
    n <- x$num_obs
    m <- matrix(0, nrow = n, ncol = n)

    for (id in 1:n) {
      nn <- x$GetNeighborSize(id - 1)
      nbrs <- x$GetNeighbors(id - 1)
      nbr_weights <- x$GetNeighborWeights(id - 1)

      for (i in 1:nn) {
        nid <- nbrs[i] + 1
        wv <- nbr_weights[i]
        m[id, nid] <- wv
      }
    }

    return(m)
  }
}

#################################################################
#' @title Read a .GAL file
#' @description Create a spatial weights object from a .GAL file
#' @param file_path The file paht of the .GAL file
#' @param id_vec The id_vec is the id values used in the .GAL file.
#' Default is empty.
#' @return A weights object
#' @export
read_gal <- function(file_path, id_vec = c()) {

  # read first line from file
  con <- file(file_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)

  items <- strsplit(first_line, " ")
  num_obs <- as.numeric(items[[1]][[2]])

  if (length(id_vec) == 0) {
    id_vec <- 0 : num_obs - 1
  }

  if (class(id_vec) == "numeric") {
    id_vec <- as.character(id_vec)
  }

  w <- p_gda_load_gal(file_path, id_vec)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title Read a .GWT file
#' @description Create a spatial weights object from a .GWT file
#' @param file_path The file paht of the .GWT file
#' @param id_vec The id_vec is the id values used in the .GWT file.
#' Default is empty.
#' @return A weights object
#' @export
read_gwt <- function(file_path, id_vec = c()) {
  # read first line from file
  con <- file(file_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)

  items <- strsplit(first_line, " ")
  num_obs <- as.numeric(items[[1]][[2]])

  if (length(id_vec) == 0) {
    id_vec <- 0 : num_obs - 1
  }

  if (class(id_vec) == "numeric") {
    id_vec <- as.character(id_vec)
  }

  w <- p_gda_load_gwt(file_path, id_vec)

  return(Weight$new(p_GeoDaWeight(w)))
}

#################################################################
#' @title Read a .SWM file
#' @description Create a spatial weights object from a .SWM file
#' @param file_path The file paht of the .SWM file
#' @param id_vec The id_vec is the id values used in the .SWM file.
#'  e.g. c(0,1,2,3,...)
#' @return A weights object
#' @export
read_swm <- function(file_path, id_vec = numeric()) {
  w <- p_gda_load_swm(file_path, id_vec)

  return(Weight$new(p_GeoDaWeight(w)))
}
