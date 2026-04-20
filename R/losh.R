#################################################################
#' @title Local Spatial Heteroscedasticity (LOSH)
#' @description The function to apply local spatial heteroscedasticity (LOSH)
#' statistics. It identifies areas of local spatial heteroscedasticity,
#' which can be used to characterize the internal heterogeneity of hotspots.
#' The cluster indicators are: 0 (Not significant), 1 (Heterogeneous), 
#' and 2 (Homogeneous). Indices 5 (Undefined) and 6 (Isolated) are also used.
#' @param w An instance of Weight object
#' @param df A data frame with only selected variable. E.g. guerry["Crm_prs"]
#' @param a (optional) The power (or exponent) of the absolute difference.
#' Default is 2.
#' @param permutations (optional) The number of permutations for the LISA
#' computation
#' @param permutation_method (optional) The permutation method used for the
#' LISA computation. Options are ('complete', 'lookup'). Default is 'complete'.
#' @param significance_cutoff  (optional) A cutoff value for significance
#' p-values to filter not-significant clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel
#' LISA computation
#' @param seed (optional) The seed for random number generator
#' @return An instance of LISA-class
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' losh <- local_losh(queen_w, guerry["Crm_prs"])
#' Hi <- lisa_values(losh)
#' }
#' @export
local_losh <- function(w, df, a = 2, permutations=999, permutation_method="complete",
                       significance_cutoff=0.05, cpu_threads=6,
                       seed=123456789) {
  if (w$num_obs <= 0) {
    stop("Weights object is not valid.")
  }

  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  data <- df[[1]]
  
  if (sum(is.na(data)) == length(data)) {
    stop("The input data is all NA.")
  }

  lisa_obj <- p_locallosh(w$GetPointer(), data, permutations,
                          permutation_method, significance_cutoff,
                          cpu_threads, seed, a)
                          
  return(LISA$new(p_LISA(lisa_obj)))
}
