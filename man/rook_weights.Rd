% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{rook_weights}
\alias{rook_weights}
\title{Rook Contiguity Spatial Weights}
\usage{
rook_weights(
  sf_obj,
  order = 1,
  include_lower_order = FALSE,
  precision_threshold = 0
)
}
\arguments{
\item{sf_obj}{An sf (simple feature) object}

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
library(sf)
guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry <- st_read(guerry_path)
rook_w <- rook_weights(guerry)
summary(rook_w)
}
