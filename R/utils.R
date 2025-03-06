############################################################
#' @title Natural Breaks (Jenks)
#' @description Natural Breaks group data whose boundaries are set where there
#' are relatively big differences.
#' @param k A numeric value indicates how many breaks
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' natural_breaks(k=5, guerry['Crm_prs'])
#' @export
natural_breaks <- function(k, df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_naturalbreaks(k, data))
}

############################################################
#' @title Quantile Breaks
#' @description Quantile breaks data into groups that each have the same number
#' of observations
#' @param k A numeric value indicates how many breaks
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' quantile_breaks(k=5, guerry['Crm_prs'])
#' @export
quantile_breaks <- function(k, df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_quantilebreaks(k, data))
}

############################################################
#' @title (Box) Hinge15 Breaks
#' @description Hinge15 breaks data into 6 groups like box plot groups
#' (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=1.5
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' hinge15_breaks(guerry['Crm_prs'])
#' @export
hinge15_breaks <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_hinge15breaks(data))
}

############################################################
#' @title (Box) Hinge30 Breaks
#' @description Hinge30 breaks data into 6 groups like box plot groups
#' (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=3.0
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' hinge30_breaks(guerry['Crm_prs'])
#' @export
hinge30_breaks <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_hinge30breaks(data))
}

############################################################
#' @title Percentile Breaks
#' @description Percentile breaks data into 6 groups: the lowest 1%, 1-10%,
#' 10-50%, 50-90%, 90-99% and the top 1%.
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' percentile_breaks(guerry['Crm_prs'])
#' @export
percentile_breaks <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_percentilebreaks(data))
}

############################################################
#' @title Standard Deviation Breaks
#' @description Standard deviation breaks first transforms data to standard
#' deviation units (mean=0, stddev=1), and then divide the range of values into
#' 6 groups.
#' @param df A data frame with selected variable. E.g. guerry["Crm_prs"]
#' @return A vector of numeric values of computed breaks
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' stddev_breaks(guerry['Crm_prs'])
#' @export
stddev_breaks <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  data <- df[[1]]
  return(p_stddevbreaks(data))
}

#################################################################
#' @title  Empirical Bayes(EB) Rate
#' @description The function to compute EB Rate from an event variable and a
#' base variable.
#' @param df A data frame with two selected variable: one is "event", anothor is
#'  "base" variable. E.g. guerry[c("hr60", "po60")]
#' @return A data.frame with two columns "EB Rate" and "IsNull".
#' @examples
#' \dontrun{
#' library(sf)
#' nat <- st_read("natregimes.shp")
#' ebr <- eb_rates(nat[c("HR60", "PO60")])
#' ebr
#' }
#' @export
eb_rates <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  event_data <- df[[1]]
  base_data <- df[[2]]
  return(p_eb_rate(event_data, base_data))
}

#################################################################
#' @title  Empirical Bayes(EB) Rate Standardization
#' @description The function to compute EB Rate Standardization from 
#' an event variable and a base variable.
#' @param df A data frame with two selected variable: one is "event", anothor is
#'  "base" variable. E.g. guerry[c("hr60", "po60")]
#' @return A data.frame with two columns "EB Rate" and "IsNull".
#' @examples
#' \dontrun{
#' library(sf)
#' nat <- st_read("natregimes.shp")
#' ebr <- eb_rates_standardization(nat[c("HR60", "PO60")])
#' ebr
#' }
#' @export
eb_rates_standardization <- function(df) {
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }
  event_data <- df[[1]]
  base_data <- df[[2]]
  return(p_eb_rate_standardization(event_data, base_data))
}