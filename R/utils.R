############################################################
#' @title Demean standardize
#' @description The mean for each variable is subtracting from each observation resulting in mean zero.
#' @param data An input data for median absolute deviation
#' @return A data array
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' crm_prp <- guerry_df['Crm_prp'][,1] # get values of variable "crm_prp"
#' data <- list(crm_prs, crm_prp)
#' demean(data)
#' @export
demean <- function(data) {
  return (p_gda_demean(data))
}

############################################################
#' @title Standardize Data (Z)
#' @description  Standarize data by transforming data to have zero mean and unit variance
#' @param data An input data for median absolute deviation
#' @return A data array
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' crm_prp <- guerry_df['Crm_prp'][,1] # get values of variable "crm_prp"
#' data <- list(crm_prs, crm_prp)
#' standardize(data)
#' @export
standardize <- function(data) {
  return (p_gda_standardize(data))
}

############################################################
#' @title Median absolute deviation
#' @description Median absolute deviation to measure  measure of the variability of a univariate sample of quantitative data
#' @param data An input data for median absolute deviation
#' @return A data array
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' crm_prp <- guerry_df['Crm_prp'][,1] # get values of variable "crm_prp"
#' data <- list(crm_prs, crm_prp)
#' median_absolute_deviation(data)
#' @export
median_absolute_deviation <- function(data) {
  return (p_gda_mad(data))
}

############################################################
#' @title Natural Breaks (Jenks)
#' @description Natural Breaks group data whose boundaries are set where there are relatively big differences.
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' natural_breaks(k=5, data=crm_prs)
#' @export
natural_breaks <- function(k, data) {
  return (p_naturalbreaks(k, data))
}

############################################################
#' @title Quantile Breaks
#' @description Quantile breaks data into groups that each have the same number of observations
#' @param k A numeric value indicates how many breaks
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' quantile_breaks(k=5, data=crm_prs)
#' @export
quantile_breaks <- function(k, data) {
  return (p_quantilebreaks(k, data))
}

############################################################
#' @title (Box) Hinge15 Breaks
#' @description Hinge15 breaks data into 6 groups like box plot groups (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=1.5
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' hinge15_breaks(data=crm_prs)
#' @export
hinge15_breaks <- function(data) {
  return (p_hinge15breaks(data))
}

############################################################
#' @title (Box) Hinge30 Breaks
#' @description Hinge30 breaks data into 6 groups like box plot groups (Lower outlier, < 25%, 25-50%, 50-75%, >75%, Upper outlier) with hinge=3.0
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' hinge30_breaks(data=crm_prs)
#' @export
hinge30_breaks <- function(data) {
  return (p_hinge30breaks(data))
}

############################################################
#' @title Percentile Breaks
#' @description Percentile breaks data into 6 groups: the lowest 1%, 1-10%, 10-50%, 50-90%, 90-99% and the top 1%.
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' percentile_breaks(data=crm_prs)
#' @export
percentile_breaks <- function(data) {
  return (p_percentilebreaks(data))
}

############################################################
#' @title Standard Deviation Breaks
#' @description Standard deviation breaks first transforms data to standard deviation units (mean=0, stddev=1), and then divide the range of values into 6 groups.
#' @param data An input data
#' @return A array contains values of computed breaks
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' crm_prs <- guerry_df['Crm_prs'][,1] # get values of variable "crm_prs"
#' stddev_breaks(data=crm_prs)
#' @export
stddev_breaks <- function(data) {
  return (p_stddevbreaks(data))
}
