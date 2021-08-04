#' @importFrom Rcpp evalCpp
#' @importFrom methods is new
#' @importFrom digest digest
#' @useDynLib rgeoda
NULL

.onUnload <- function(libname, pkgname) {
}

# When creating weights from sf/sp object, the geoda object will be created and
# cached. The geoda object will not be visible to users. The following cacheEnv
# is defined for internally use only.
# (1) assign(key, value, envir=cacheEnv)
# (2) get(key, envir=cacheEnv)
# key is a md5 string created using the sf/sp object. E.g.
# > digest::digest(pol_pres15)
# value is the geoda object
cacheEnv <- new.env()

# This function is used internally to get a geoda object from a sf object
# The geoda object will be created internally and will not be visible.
getGeoDaObj <- function(sf_obj) {
  # check if there is a geoda object being created in the cacheEnv
  geom_col <- sf_obj[[attr(sf_obj, "sf_column")]]
  gda_uid <- digest::digest(geom_col) # use geometry to create uid

  geoda_obj <- NULL

  if (exists(gda_uid, envir = cacheEnv)) {
    # use it directly for weights creation
    geoda_obj <- get(gda_uid, envir = cacheEnv)
  } else {
    # create a geoda object from sf_obj, NOT copy table
    geoda_obj <- as.geoda(sf_obj, with_table = FALSE)
    # cache it
    assign(gda_uid, geoda_obj, envir = cacheEnv)
  }

  if (is.null(geoda_obj) || geoda_obj$GetNumObs() <= 0) {
    stop("Failed to read sf object. Please submit a bug report to https://github.com/geodacenter/rgeoda")
  }

  return(geoda_obj)
}
