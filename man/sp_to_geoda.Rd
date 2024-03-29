% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf_geoda.R
\name{sp_to_geoda}
\alias{sp_to_geoda}
\title{Create an instance of geoda-class from a 'sp' object}
\usage{
sp_to_geoda(sp_obj, with_table = TRUE)
}
\arguments{
\item{sp_obj}{An instance of 'sp' object}

\item{with_table}{A boolean flag indicates if table is copied from sf object
to create geoda object. Default is TRUE}
}
\value{
An instance of geoda-class
}
\description{
Create an instance of geoda-class from a 'sp' object. NOTE: The
table content is NOT used to create an instance of geoda-class.
}
