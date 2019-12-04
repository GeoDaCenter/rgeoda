#################################################################
#' @title LISA class
#' @description A LISA class wrappers the results of LISA computation
#' @field gda_lisa An object of GeoDaLISA
#' @field p_vals The pseudo-p values of significance of LISA computation
#' @field c_vals The cluster indicators of LISA computation
#' @field lisa_vals The local spatial autocorrelation values of LISA computation
#' @field nn_vals The number of neighbors of every observations in LISA computation
#' @field lag_vals The values of spatial lags of LISA computation
#' @field num_obs Number of observations
#' @field labels The cluster labels of LISA
#' @field colors The cluster colors (HEX format) of LISA
#' @export
LISA <- setRefClass("LISA",
  fields = list(
    gda_lisa = "_p_LISA",
    p_vals = "array",
    c_vals = "array",
    lisa_vals = "array",
    nn_vals = "array",
    lag_vals = "array",
    num_obs = "integer",
    labels = "array",
    colors = "array"
  ),
  methods = list(
    initialize = function(lisa_obj) {
      "Constructor with a LISA object (internally used)"
      .self$gda_lisa = lisa_obj
    },
    Run = function() {
      "Call to run LISA computation"
      gda_lisa$Run()
    },
    SetPermutations = function(num_perm) {
      "Set the number of permutations for the LISA computation"
      if (num_perm < 1 || num_perm > 999999) {
        stop("The number of permutations is a positive integer number, but has to be less than 999999.")
      }
      gda_lisa$SetNumPermutations(num_perm)
    },
    SetThreads= function(num_threads) {
      "Set the number of CPU threads for the LISA computation"
      if (num_threads < 1) {
        stop("The number of CPU threads has to be a positive integer number.")
      }
      gda_lisa$SetNumThreads(num_threads)
    },
    GetLISAValues = function() {
      "Get the local spatial autocorrelation values returned from LISA computation."
      return (gda_lisa$GetLISAValues())
    },
    GetPValues = function() {
      "Get the local pseudo-p values of significance returned from LISA computation."
      return (gda_lisa$GetLocalSignificanceValues())
    },
    GetClusterIndicators = function() {
      "Get the local cluster indicators returned from LISA computation."
      return (gda_lisa$GetClusterIndicators())
    },
    GetNumNeighbors = function() {
      "Get the number of neighbors of every observations in LISA computation."
      return (gda_lisa$GetNumNeighbors())
    },
    GetLabels = function() {
      "Get the cluster labels of LISA computation."
      return (gda_lisa$GetLabels())
    },
    GetColors = function() {
      "Get the cluster colors of LISA computation."
      return (gda_lisa$GetColors())
    }
  )
)

#################################################################
#' @title  Local Moran statistics
#' @description The function to apply local Moran statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_moran <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_localmoran(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Geary statistics
#' @description The function to apply local Geary statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_geary <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_geary(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Multivariate Geary statistics
#' @description The function to apply local Multivariate Geary statistics
#' @param w An instance of Weight object
#' @param data A 2D tuple of values of selected variables
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_multigeary <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) <  1) {
    stop("The number of variables has to be larger than 1.")
  }

  lisa_obj <- gda_multigeary(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Getis-Ord's G statistics
#' @description The function to apply Getis-Ord's local G statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_g <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_localg(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Getis-Ord's G* statistics
#' @description The function to apply Getis-Ord's local G* statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_gstar <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_localgstar(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Join Count statistics
#' @description The function to apply local Join Count statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_joincount <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_joincount(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Local Multivariate Join Count statistics
#' @description The function to apply local Multivariate Join Count statistics
#' @param w An instance of Weight object
#' @param data A 2D tuple of values of selected variables
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_multijoincount <- function(w, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) <  1) {
    stop("The number of variables has to be larger than 1.")
  }

  lisa_obj <- gda_multijoincount(w$gda_w, data)
  return (LISA$new(lisa_obj))
}

#################################################################
#' @title  Quantile LISA statistics
#' @description The function to apply quantile LISA statistics
#' @param w An instance of Weight object
#' @param data A tuple of values of selected variable
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @export
local_quantileSA <- function(w, k, data) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- gda_quantilelisa(w$gda_w, k, data)
  return (LISA$new(lisa_obj))
}
