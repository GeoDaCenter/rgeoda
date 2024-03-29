% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{gda_rook_weights}
\alias{gda_rook_weights}
\title{(For internally use and test only) Rook Contiguity Spatial Weights}
\usage{
gda_rook_weights(
  geoda_obj,
  order = 1,
  include_lower_order = FALSE,
  precision_threshold = 0
)
}
\arguments{
\item{geoda_obj}{An object of [geoda] class}

\item{order}{(Optional) Order of contiguity}

\item{include_lower_order}{(Optional)  Whether or not the lower order
neighbors should be included in the weights structure}

\item{precision_threshold}{(Optional) The precision of the underlying shape
file is insufficient to allow for an exact match of coordinates to determine
 which polygons are neighbors}
}
\value{
An instance of Weight-class
}
\description{
Create a Rook contiguity weights with options of "order",
"include lower order" and "precision threshold"
}
\examples{
\dontrun{
guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry <- geoda_open(guerry_path)
rook_w <- gda_rook_weights(guerry)
summary(rook_w)
}
}
