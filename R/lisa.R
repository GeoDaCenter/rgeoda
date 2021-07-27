#################################################################
#' @title LISA class (Internally Used)
#' @description A LISA-class that wrappers the statistics of LISA computation
#' @field gda_lisa An object of GeoDaLISA
#' @field p_vals The pseudo-p values of significance of LISA computation
#' @field c_vals The cluster indicators of LISA computation
#' @field lisa_vals The local spatial autocorrelation values of LISA computation
#' @field nn_vals The number of neighbors of every observations in LISA computation
#' @field labels The cluster labels of LISA
#' @field colors The cluster colors (HEX format) of LISA
#' @export
LISA <- setRefClass("LISA",
  fields = list(
    gda_lisa = "p_LISA",
    p_vals = "numeric",
    c_vals = "numeric",
    lisa_vals = "numeric",
    nn_vals = "numeric",
    labels = "character",
    colors = "character"
  ),
  methods = list(
    initialize = function(lisa_obj) {
      "Constructor with a LISA object (internally used)"
      .self$gda_lisa = lisa_obj
      .self$p_vals = .self$GetLocalSignificanceValues()
      .self$c_vals = lisa_obj$GetClusterIndicators()
      .self$lisa_vals = lisa_obj$GetLISAValues()
      .self$nn_vals = lisa_obj$GetNumNeighbors()
      .self$labels = lisa_obj$GetLabels()
      .self$colors = lisa_obj$GetColors()
    },
    Run = function() {
      "Call to run LISA computation"
      gda_lisa$Run()
      # update values
      .self$p_vals = .self$GetLocalSignificanceValues()
      .self$c_vals = gda_lisa$GetClusterIndicators()
      .self$lisa_vals = gda_lisa$GetLISAValues()
      .self$nn_vals = gda_lisa$GetNumNeighbors()
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
      pvals <- gda_lisa$GetLocalSignificanceValues()
      num_obs <- length(pvals)
      for (row_idx in 1:num_obs) {
        if (pvals[row_idx] == -1) {
          pvals[row_idx] <- NA
        }
      }
      return (pvals)
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
#' @title  Bonferroni bound value of local spatial autocorrelation
#' @description Get Bonferroni bound value based on current LISA computation and current significat p-value
#' @param gda_lisa An instance of LISA object
#' @param current_p A value of current siginificant p-value
#' @return A numeric value of Bonferroni bound
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' bo <- lisa_bo(lisa, 0.05)
#' bo
#' }
#' @export
lisa_bo <- function(gda_lisa, current_p) {
  return (gda_lisa$GetBO(current_p))
}

#################################################################
#' @title  False Discovery Rate value of local spatial autocorrelation
#' @description Get False Discovery Rate value based on current LISA computation and current significant p-value
#' @param gda_lisa An instance of LISA object
#' @param current_p A value of current siginificant p-value
#' @return A numeric vector of False Discovery Rate
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' fdr <- lisa_fdr(lisa, 0.05)
#' fdr
#' }
#' @export
lisa_fdr <- function(gda_lisa, current_p) {
  return (gda_lisa$GetFDR(current_p))
}

#################################################################
#' @title  Get LISA values
#' @description Get the local spatial autocorrelation values returned from LISA computation
#' @param gda_lisa An instance of LISA object
#' @return A numeric vector of local spatial autocorrelation
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' lms <- lisa_values(lisa)
#' lms
#' }
#' @export
lisa_values <- function(gda_lisa) {
  return (gda_lisa$GetLISAValues())
}

#################################################################
#' @title  Get pseudo-p values of LISA
#' @description Get the local pseudo-p values of significance returned from LISA computation.
#' @param gda_lisa An instance of LISA object
#' @return A numeric vector of pseudo-p values of local spatial autocorrelation
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' pvals <- lisa_pvalues(lisa)
#' pvals
#' }
#' @export
lisa_pvalues <- function(gda_lisa) {
  return (gda_lisa$GetLocalSignificanceValues())
}

#################################################################
#' @title  Get local cluster indicators
#' @description Get the local cluster indicators returned from LISA computation.
#' @param gda_lisa An instance of LISA object
#' @param cutoff A value of cutoff for significance p-values to filter not-significant clusters, default=0.0, means not used
#' @return A numeric vector of LISA cluster indicator
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' clsts <- lisa_clusters(lisa)
#' clsts
#' }
#' @export
lisa_clusters <- function(gda_lisa, cutoff=0.0) {
  if (cutoff > 0.0) {
    gda_lisa$SetSignificanceCutoff(cutoff)
  }
  return (gda_lisa$GetClusterIndicators())
}

#################################################################
#' @title  Get numbers of neighbors for all observations
#' @description Get numbers of neighbors for all observations
#' @param gda_lisa An instance of LISA object
#' @return A numeric vector of the number of neighbors
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
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
#' @return A string vector of cluster labels
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
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
#' @return A string vector of cluster colors
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' clrs <- lisa_colors(lisa)
#' clrs
#' }
lisa_colors <- function(gda_lisa) {
  return (gda_lisa$GetColors())
}

#################################################################
#' @title  Local Moran Statistics
#' @description The function to apply local Moran statistics
#' @param w An instance of Weight object
#' @param df A data frame with only selected variable. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_moran(queen_w, guerry["Crm_prs"])
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_moran <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]
  lisa_obj <- p_localmoran(w$GetPointer(), data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Moran with Empirical Bayes(EB) Rate
#' @description The function to apply local Moran with EB Rate statistics. The EB rate is first computed from "event" and "base" variables, and then used in local moran statistics.
#' @param w An instance of Weight object
#' @param df A data frame with two selected variable: one is "event", anothor is "base" variable. E.g. guerry[c("hr60", "po60")]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' \dontrun{
#' library(sf)
#' nat <- st_read("natregimes.shp")
#' nat_w <- queen_weights(nat)
#' lisa <- local_moran_eb(queen_w, guerry[c("hr60", "po60")])
#' lms <- lisa_values(lisa)
#' lms
#' }
#' @export
local_moran_eb <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (class(w)[[1]] != "Weight") {
    stop("The parameter 'w' needs to be an instance of Weight object.")
  }
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  event_data <- df[[1]]
  base_data <- df[[2]]

  lisa_obj <- p_localmoran_eb(w$GetPointer(), event_data, base_data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Geary Statistics
#' @description The function to apply local Geary statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_geary(queen_w, guerry["Crm_prs"])
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_geary <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]
  lisa_obj <- p_localgeary(w$GetPointer(), data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)

  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Multivariate Geary Statistics
#' @description The function to apply local Multivariate Geary statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variables only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' lisa <- local_multigeary(queen_w, data)
#' lms <- lisa_clusters(lisa)
#' lms
#' @export
local_multigeary <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  num_vars <- length(df)

  if (inherits(df, "sf")) {
    num_vars <- num_vars - 1
  }

  lisa_obj <- p_localmultigeary(w$GetPointer(), df, num_vars, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Getis-Ord's G Statistics
#' @description The function to apply Getis-Ord's local G statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_g(queen_w, guerry["Crm_prs"])
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_g <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]

  lisa_obj <- p_localg(w$GetPointer(), data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Getis-Ord's G* Statistics
#' @description The function to apply Getis-Ord's local G* statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_gstar(queen_w,  guerry["Crm_prs"])
#' lms <- lisa_values(lisa)
#' lms
#' @export
local_gstar <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]

  lisa_obj <- p_localgstar(w$GetPointer(), data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Local Join Count Statistics
#' @description The function to apply local Join Count statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_joincount(queen_w, guerry['TopCrm'])
#' clsts<- lisa_clusters(lisa)
#' clsts
#' @export
local_joincount <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]

  lisa_obj <- p_localjoincount(w$GetPointer(), data, permutations, permutation_method, significance_cutoff, cpu_threads, seed)

  jc <- LISA$new(p_LISA(lisa_obj))

  # update the probability results: change these with jc=0 to NA
  for (idx in 1:w$num_obs) {
    if (jc$lisa_vals[idx] == 0) {
      jc$p_vals[idx] <- NA
    }
  }

  return(jc)
}

#################################################################
#' @title  Bivariate Local Join Count Statistics
#' @description The function to apply local Bivariate Join Count statistics
#' @param w An instance of Weight object
#' @param df A data frame with two selected variable. E.g. guerry[c("TopCrm", "InvCrm")]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry["InvCrm"] <-  1 - guerry[["TopCrm"]]
#' lisa <- local_bijoincount(queen_w, guerry[c("TopCrm", "InvCrm")])
#' clsts<- lisa_clusters(lisa)
#' clsts
#' @export
local_bijoincount <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data1 <- df[[1]]
  data2 <- df[[2]]

  if (p_gda_isbinary(data1) == FALSE || p_gda_isbinary(data2) == FALSE) {
    stop("The input data is not binary.")
  }

  if (sum(data1 + data2) != w$num_obs) {
    stop("The bivariate local join count only applies on two variables with no-colocation.")
  }

  lisa_obj <- p_localmultijoincount(w$GetPointer(), df, 2, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  jc <- LISA$new(p_LISA(lisa_obj))

  # update the probability results: change these with jc=0 to NA
  for (idx in 1:w$num_obs) {
    if (jc$lisa_vals[idx] == 0) {
      jc$p_vals[idx] <- NA
    }
  }

  return(jc)
}

#################################################################
#' @title (Multivariate) Colocation Local Join Count Statistics
#' @description The function to apply (multivariate) colocation local Join Count statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variables only. E.g. guerry[c("TopCrm", "TopWealth", "TopLit")]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_multijoincount(queen_w,  guerry[c('TopWealth','TopWealth', 'TopLit')])
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_multijoincount <- function(w, df, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  num_vars <- length(df)

  if (inherits(df, "sf")) {
    num_vars <- num_vars - 1
  }

  for ( idx in 1:num_vars)  {
    if (p_gda_isbinary(df[[idx]]) == FALSE) {
      stop("The input data is not binary.")
    }
  }

  if (num_vars == 2) {
    if (sum(df[[1]] + df[[2]]) == w$num_obs) {
      stop("The input two variables have no colocations. Please use bivariate local join count: local_bijoincount().")
    }
  }

  lisa_obj <- p_localmultijoincount(w$GetPointer(), df, num_vars, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  jc <- LISA$new(p_LISA(lisa_obj))

  # update the probability results: change these with jc=0 to NA
  for (idx in 1:w$num_obs) {
    if (jc$lisa_vals[idx] == 0) {
      jc$p_vals[idx] <- NA
    }
  }

  return(jc)
}

#################################################################
#' @title  Quantile LISA Statistics
#' @description The function to apply quantile LISA statistics
#' @param w An instance of Weight object
#' @param k A value indicates the number of quantiles. Value range e.g. [1, 10]
#' @param q A value indicates which quantile or interval used in local join count statistics. Value stars from 1.
#' @param df A data frame with selected variable only. E.g. guerry["Crm_prs"]
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_quantilelisa(queen_w, guerry["Crm_prs"], k=4, q=1)
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_quantilelisa <- function(w, df, k, q, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  num_vars <- length(df)

  if (q < 1 || q > k) {
    stop("The value of which quantile been selected should be in the range of [1, k]")
  }

  lisa_obj <- p_quantilelisa(w$GetPointer(), k, q, df[[1]], permutations, permutation_method, significance_cutoff, cpu_threads, seed)

  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title  Multivariate Quantile LISA Statistics
#' @description The function to apply multivariate quantile LISA statistics
#' @param w An instance of Weight object
#' @param df A data frame with selected variables only. E.g. guerry[c("TopCrm", "TopWealth", "TopLit")]
#' @param k A vector of "k" values indicate the number of quantiles for each variable. Value range e.g. [1, 10]
#' @param q A vector of "q" values indicate which quantile or interval for each variable used in local join count statistics. Value stars from 1.
#' @param permutations (optional) The number of permutations for the LISA computation
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are {'complete', 'lookup'}. Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' lisa <- local_multiquantilelisa(queen_w, guerry[c("Crm_prp", "Litercy")], k=c(4,4), q=c(1,1))
#' clsts <- lisa_clusters(lisa)
#' clsts
#' @export
local_multiquantilelisa <- function(w, df, k, q, permutations=999, permutation_method="complete", significance_cutoff=0.05, cpu_threads=6, seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars <= 0) {
    stop("Please specify more than one variable for multi-quantile lisa.")
  }

  if (length(k) != length(q) || length(k) != n_vars) {
    stop("Please specify 'k' and 'q' values for each variable.")
  }

  for (i in 1:n_vars) {
    ki <- k[i]
    qi <- q[i]

    if (qi < 1 || qi > ki) {
      stop("The value of which quantile been selected should be in the range of [1, k]")
    }
  }

  lisa_obj <- p_multiquantilelisa(w$GetPointer(), k, q, df, permutations, permutation_method, significance_cutoff, cpu_threads, seed)
  return (LISA$new(p_LISA(lisa_obj)))
}

#################################################################
#' @title Local Neighbor Match Test
#' @description The local neighbor match test is to assess the extent of overlap between k-nearest neighbors in geographical space and k-nearest neighbors in multi-attribute space.
#' @param df A subset of sf object with selected variables. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param k a positive integer number for k-nearest neighbors searching.
#' @param scale_method (optional) One of the scaling methods {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'} to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The type of distance metrics used to measure the distance between input data. Options are {'euclidean', 'manhattan'}. Default is 'euclidean'.
#' @param power (optional) The power (or exponent) of a number says how many times to use the number in a multiplication.
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on distance value.
#' @param is_arc (optional) FALSE (default) or TRUE, compute arc distance between two observations.
#' @param is_mile (optional) TRUE (default) or FALSE, convert distance unit from mile to km.
#' @return A data.frame with two columns "Cardinality" and "Probability".
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' nbr_test <- neighbor_match_test(data, 6)
#' nbr_test
#' @export
neighbor_match_test <- function(df, k, scale_method = "standardize", distance_method = "euclidean", power = 1.0, is_inverse = FALSE,
                                is_arc = FALSE, is_mile = TRUE) {
  if (inherits(df, "sf") == FALSE) {
    stop("The input data needs to be a sf object.")
  }
  geoda_obj <- getGeoDaObj(df) # get geoda_obj from cache or create instantly
  n_vars <- length(df) - 1 # minus geometry column

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  result <- p_neighbor_match_test(geoda_obj$GetPointer(), k, power, is_inverse, is_arc, is_mile, df, n_vars, scale_method, distance_method)

  # update the probability results: change those with -1 to NA
  for (row_idx in 1:geoda_obj$n_obs) {
    if (result[row_idx, 2] == -1) {
      result[row_idx, 2] <- NA
    }
  }

  return(result)
}
