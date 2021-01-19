#################################################################
#' @title LISA class (Internally Used)
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
    gda_lisa = "p_LISA",
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
    GetLocalSignificanceValues = function() {
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
    SetSignificanceCutoff = function(cutoff) {
      "Set the cutoff value of significance values"
      return (gda_lisa$SetSignificanceCutoff(cutoff))
    },
    GetFDR = function(current_p) {
      "Get the False Discovery Rate value"
      return (gda_lisa$GetFDR(current_p))
    },
    GetBO = function(current_p) {
      "Get the Bonferroni bound value"
      return (gda_lisa$GetBO(current_p))
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
#' @title  Get Bonferroni bound value of local spatial autocorrelation
#' @description Get Bonferroni bound value based on current LISA computation and current significat p-value
#' @param gda_lisa An instance of LISA object
#' @param current_p A value of current siginificant p-value
#' @return A value of Bonferroni bound
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' bo <- lisa_bo(lisa, 0.05)
#' bo
#' }
#' @export
lisa_bo <- function(gda_lisa, current_p) {
  return (gda_lisa$GetBO(current_p))
}

#################################################################
#' @title  Get False Discovery Rate value of local spatial autocorrelation
#' @description Get False Discovery Rate value based on current LISA computation and current significat p-value
#' @param gda_lisa An instance of LISA object
#' @param current_p A value of current siginificant p-value
#' @return A value of False Discovery Rate
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' fdr <- lisa_fdr(lisa, 0.05)
#' fdr
#' }
#' @export
lisa_fdr <- function(gda_lisa, current_p) {
  return (gda_lisa$GetFDR(current_p))
}

#################################################################
#' @title  Get values of local spatial autocorrelation
#' @description Get the local spatial autocorrelation values returned from LISA computation
#' @param gda_lisa An instance of LISA object
#' @return data A numeric array of local spatial autocorrelation
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' lms <- lisa_values(lisa)
#' lms
#' }
#' @export
lisa_values <- function(gda_lisa) {
  return (gda_lisa$GetLISAValues())
}

#################################################################
#' @title  Get values of local pseudo-p values
#' @description Get the local pseudo-p values of significance returned from LISA computation.
#' @param gda_lisa An instance of LISA object
#' @return data a tuple of pseudo-p values of local spatial autocorrelation
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' pvals <- lisa_pvalues(lisa)
#' pvals
#' }
#' @export
lisa_pvalues <- function(gda_lisa) {
  return (gda_lisa$GetLocalSignificanceValues())
}

#################################################################
#' @title  Get values of local cluster indicators
#' @description Get the local cluster indicators returned from LISA computation.
#' @param gda_lisa An instance of LISA object
#' @param cutoff A value of cutoff for significance p-values to filter not-significant clusters, default=0.05
#' @return data a tuple values
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' clsts <- lisa_clusters(lisa)
#' clsts
#' }
#' @export
lisa_clusters <- function(gda_lisa, cutoff=0.05) {
  gda_lisa$SetSignificanceCutoff(cutoff)
  return (gda_lisa$GetClusterIndicators())
}

#################################################################
#' @title  Get numbers of neighbors for all observations
#' @description Get numbers of neighbors for all observations
#' @param gda_lisa An instance of LISA object
#' @return data A numeric array
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' nn <- lisa_num_nbrs(lisa)
#' nn
#' }
#' @export
lisa_num_nbrs <- function(gda_lisa) {
  return (gda_lisa$GetNumNeighbors())
}

#################################################################
#' @title  Get cluster labels
#' @description Get cluster labels of LISA computation.
#' @param gda_lisa An instance of LISA object
#' @return labels A numeric array
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' lbls <- lisa_labels(lisa)
#' lbls
#' }
#' @export
lisa_labels <- function(gda_lisa) {
  return (gda_lisa$GetLabels())
}

#################################################################
#' @title  Get cluster colors
#' @description Get the cluster colors of LISA computation.
#' @param gda_lisa An instance of LISA object
#' @return colors A numeric array
#' @export
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' clrs <- lisa_colors(lisa)
#' clrs
#' }
lisa_colors <- function(gda_lisa) {
  return (gda_lisa$GetColors())
}

#################################################################
#' @title  Local Moran statistics
#' @description The function to apply local Moran statistics
#' @param w An instance of Weight object
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_moran(queen_w, crm_prs)
#' lms <- lisa_values(gda_lisa = lisa)
#' lms
#' @export
local_moran <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- p_localmoran(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Geary statistics
#' @description The function to apply local Geary statistics
#' @param w An instance of Weight object
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_geary(queen_w, crm_prs)
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_geary <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- p_localgeary(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Multivariate Geary statistics
#' @description The function to apply local Multivariate Geary statistics
#' @param w An instance of Weight object
#' @param data A list of numeric values of selected variables
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' lisa <- local_multigeary(queen_w, data)
#' lms <- lisa_clusters(lisa)
#' lms
#' @export
local_multigeary <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) <  1) {
    stop("The number of variables has to be larger than 1.")
  }

  lisa_obj <- p_localmultigeary(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Getis-Ord's G statistics
#' @description The function to apply Getis-Ord's local G statistics
#' @param w An instance of Weight object
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_g(queen_w, crm_prs)
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_g <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- p_localg(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Getis-Ord's G* statistics
#' @description The function to apply Getis-Ord's local G* statistics
#' @param w An instance of Weight object
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' lisa <- local_gstar(queen_w, crm_prs)
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_gstar <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- p_localgstar(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Join Count statistics
#' @description The function to apply local Join Count statistics
#' @param w An instance of Weight object
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' top_crm <- guerry_df['TopCrm'][,1] # get 0/1 values of variable "top_crm"
#' lisa <- local_joincount(queen_w, top_crm)
#' clsts<- lisa_clusters(lisa)
#' clsts
#' @export
local_joincount <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }

  lisa_obj <- p_localjoincount(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Bivariate Join Count statistics
#' @description The function to apply local Bivariate Join Count statistics
#' @param w An instance of Weight object
#' @param data1 A numeric array of selected variable
#' @param data2 A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' top_crm <- guerry_df['TopCrm'][,1]
#' inv_crm <-  1 - top_crm
#' lisa <- local_bijoincount(queen_w, top_crm, inv_crm)
#' clsts<- lisa_clusters(lisa)
#' clsts
#' @export
local_bijoincount <- function(w, data1, data2, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data1) != w$num_obs) {
    stop("The size of data1 doesnt not match the number of observations")
  }
  if (length(data2) != w$num_obs) {
    stop("The size of data2 doesnt not match the number of observations")
  }
  if (sum(data1 + data2) != w$num_obs) {
    stop("The bivariate local join count only applies on two variables with no-colocation.")
  }
  data <- list(data1, data2)
  lisa_obj <- p_localmultijoincount(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Multivariate Join Count statistics
#' @description The function to apply local Multivariate Join Count statistics
#' @param w An instance of Weight object
#' @param data A list of numeric values of selected variables
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' bin_data <- guerry_df[c('TopWealth','TopWealth', 'TopLit')]
#' lisa <- local_multijoincount(queen_w, bin_data)
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_multijoincount <- function(w, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) <  1) {
    stop("The number of variables has to be larger than 1.")
  }

  lisa_obj <- p_localmultijoincount(w$GetPointer(), data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Quantile LISA statistics
#' @description The function to apply quantile LISA statistics
#' @param w An instance of Weight object
#' @param k A value indicates the number of quantiles. Value range e.g. [1, 10]
#' @param q A value indicates which quantile or interval used in local join count statistics. Value stars from 1.
#' @param data A numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1]
#' lisa <- local_quantilelisa(queen_w, k=4, q=1, crm_prs)
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_quantilelisa <- function(w, k, q, data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (length(data) != w$num_obs) {
    stop("The size of data doesnt not match the number of observations")
  }
  if (q < 1 || q > k) {
    stop("The value of which quantile been selected should be in the range of [1, k]")
  }

  lisa_obj <- p_quantilelisa(w$GetPointer(), k, q, data, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Multivariate Quantile LISA statistics
#' @description The function to apply multivariate quantile LISA statistics
#' @param w An instance of Weight object
#' @param quantile_data A list of [k, q, data] for more than one variable. Each variable will be set with: k, indicates the number of quantiles; q, indicates which quantile or interval used in local join count statistics; data, is a numeric array of selected variable
#' @param permutations The number of permutations for the LISA computation
#' @param significance_cutoff  A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads The number of cpu threads used for parallel LISA computation
#' @param seed The seed for random number generator
#' @return lisa_obj An instance of LISA (LocalSpatialAutocorrelation) object
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prp <- guerry_df['Crm_prp'][,1]
#' lit <- guerry_df['Litercy'][,1]
#' quantiles <- list(list(4,1,crm_prp), list(4,1, lit))
#' lisa <- local_multiquantilelisa(queen_w, quantiles)
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_multiquantilelisa <- function(w, quantile_data, permutations=999, significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  n_vars = length(quantile_data)

  if (n_vars <= 0) {
    stop("Please specify quantile for at least one variable.")
  }

  k_s <- vector('numeric', n_vars)
  q_s <- vector('numeric', n_vars)
  data_s <- vector('list', n_vars)

  for (i in 1:n_vars) {
    qd <- quantile_data[[i]]
    k <- qd[[1]]
    q <- qd[[2]]

    if (q < 1 || q > k) {
      stop("The value of which quantile been selected should be in the range of [1, k]")
    }

    data <- qd[[3]]
    k_s[[i]] <- k
    q_s[[i]] <- q
    data_s[[i]] <- data
  }

  lisa_obj <- p_multiquantilelisa(w$GetPointer(), k_s, q_s, data_s, permutations, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}
