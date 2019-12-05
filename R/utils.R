
demean <- function(data) {
  return (gda_demean(data))
}

standardize <- function(data) {
  return (gda_standardize(data))
}

mad <- function(data) {
  return (gda_standardize_mad(data))
}

natural_breaks <- function(k, data) {
  return (gda_naturalbreaks(k, data))
}

quantile_breaks <- function(k, data) {
  return (gda_guantilebreaks(k, data))
}

hinge15_breaks <- function(k, data) {
  return (gda_hinge15breaks(k, data))
}

hinge30_breaks <- function(k, data) {
  return (gda_hinge30breaks(k, data))
}

percentile_breaks <- function(k, data) {
  return (gda_percentilebreaks(k, data))
}

stddev_breaks <- function(k, data) {
  return (gda_stddevbreaks(k, data))
}
