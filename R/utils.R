
############################################################
#' @title Demean standardize
#' @description The mean for each variable is subtracting from each observation resulting in mean zero.
#' @param data An input data for median absolute deviation
#' @return A data array
#' @export
demean <- function(data) {
  return (gda_demean(data))
}

############################################################
#' @title Standardize Data (Z)
#' @description  Standarize data by transforming data to have zero mean and unit variance
#' @param data An input data for median absolute deviation
#' @return A data array
#' @export
standardize <- function(data) {
  return (gda_standardize(data))
}

############################################################
#' @title Median absolute deviation
#' @description Median absolute deviation to measure  measure of the variability of a univariate sample of quantitative data
#' @param data An input data for median absolute deviation
#' @return A data array
#' @export
median_absolute_deviation <- function(data) {
  return (gda_standardize_mad(data))
}

############################################################
#' @title Natural Breaks (Jenks)
#' @description Natural Breaks group data whose boundaries are set where there are relatively big differences.
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
natural_breaks <- function(k, data) {
  return (gda_naturalbreaks(k, data))
}

############################################################
#' @title Quantile Breaks
#' @description Quantile breaks data into groups that each have the same number of observations
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
quantile_breaks <- function(k, data) {
  return (gda_guantilebreaks(k, data))
}

############################################################
#' @title (Box) Hinge15 Breaks
#' @description Hinge15 breaks data into groups like box plot groups (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=1.5
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
hinge15_breaks <- function(k, data) {
  return (gda_hinge15breaks(k, data))
}

############################################################
#' @title (Box) Hinge30 Breaks
#' @description Hinge30 breaks data into groups like box plot groups (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=3.0
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
hinge30_breaks <- function(k, data) {
  return (gda_hinge30breaks(k, data))
}

############################################################
#' @title Percentile Breaks
#' @description Percentile breaks data into 6 groups: the lowest 1%, 1-10%, 10-50%, 50-90%, 90-99% and the top 1%.
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
percentile_breaks <- function(k, data) {
  return (gda_percentilebreaks(k, data))
}

############################################################
#' @title Standard Deviation Breaks
#' @description Standard deviation breaks first transforms data to standard deviation units (mean=0, stddev=1), and then divide the range of values into 6 groups.
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @export
stddev_breaks <- function(k, data) {
  return (gda_stddevbreaks(k, data))
}
