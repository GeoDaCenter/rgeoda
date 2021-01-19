# This file is used to wrap C++ classes and functions defines in RcppExports.R
# All other R script files will use this file as a bridge to C++ classes and functions
#
# Author: lixun910@gmail.com
# Changes:
# 10/29/2020 Add p_GeoDa class
# 12/23/2020 Add p_Weight class
# 12/23/2020 Add p_LISA class


#' @title p_GeoDa
#' @description p_GeoDa class is a RefClass that wraps the C++ GeoDa class
#' See C++ functions in rcpp_rgeoda.cpp
#' @export
p_GeoDa <- setClass( "p_GeoDa", representation( pointer = "externalptr" ) )

# GeoDa_method, helper function to generate C functions
# e.g. p_GeoDa__GetNumObs
# Methods are listed in RcppExports.R
#
p_GeoDa_method <- function(name) {
  paste( "p_GeoDa", name, sep = "__" )
}

#' @name $,p_GeoDa-method
#' @aliases $,p_GeoDa-method
#' @docType methods
#' @rdname p_GeoDa-class
NULL
setMethod( "$", "p_GeoDa", function(x = "p_GeoDa", name = "ANY") {
  function(...) do.call( p_GeoDa_method(name) , list(x@pointer , ... ))
})

# Constructors for p_GeoDa class
# Note: here simply using argc to determine which constructor should be called
#
setMethod( "initialize", "p_GeoDa", function(.Object, ...) {
  argtypes <- mapply(class, list(...));
  #argv <- list(...);
  argc <- length(argtypes);

  if (argc == 0) {
    # this is for uring p_GeoDa as a member in class('geoda')  in read_geoda.R
  } else if (argc > 1)  {
    # create p_GeoDa instance from sf/sp object
    .Object@pointer <- do.call( p_GeoDa_method("new1"), list(...) )
  } else {
    # create p_GeoDa instance from ESRI shapefile
    .Object@pointer <- do.call( p_GeoDa_method("new"), list(...) )
  }
  .Object
})



#' @title p_GeoDaWeight
#' @description p_GeoDaWeight class is a RefClass that wraps the C++ GeoDaWeight class
#' See C++ functions in rcpp_weights.cpp
#' @export
p_GeoDaWeight <- setClass( "p_GeoDaWeight", representation( pointer = "externalptr" ) )

# GeoDa_method, helper function to generate C functions
# e.g. p_GeoDa__GetNumObs
# Methods are listed in RcppExports.R
#
p_GeoDaWeight_method <- function(name) {
  paste( "p_GeoDaWeight", name, sep = "__" )
}

#' @name $,p_GeoDaWeight-method
#' @aliases $,p_GeoDaWeight-method
#' @docType methods
#' @rdname p_GeoDaWeight-class
NULL
setMethod( "$", "p_GeoDaWeight", function(x = "p_GeoDaWeight", name = "ANY") {
  function(...) do.call( p_GeoDaWeight_method(name) , list(x@pointer , ... ))
})

# Constructors for p_GeoDaWeight class
# Note: here simply using argc to determine which constructor should be called
#
setMethod( "initialize", "p_GeoDaWeight", function(.Object, ...) {
  argv = list(...)
  argtypes <- mapply(class, argv);
  argc <- length(argtypes);

  if (argc == 0) {
    # this is for using p_GeoDaWeight as a member in class('weight')  in weights.R
  } else {
    .Object@pointer <- argv[[1]]
  }
  .Object
})


#' @title p_LISA
#' @description p_LISA class is a RefClass that wraps the C++ LISA class
#' See C++ functions in rcpp_lisa.cpp
#' @export
p_LISA <- setClass( "p_LISA", representation( pointer = "externalptr" ) )

# LISA_method, helper function to generate C functions
# e.g. p_LISA_Run
# Methods are listed in RcppExports.R
#
p_LISA_method <- function(name) {
  paste( "p_LISA", name, sep = "__" )
}

#' @name $,p_LISA-method
#' @aliases $,p_LISA-method
#' @docType methods
#' @rdname p_LISA-class
NULL
setMethod( "$", "p_LISA", function(x = "p_LISA", name = "ANY" ) {
  function(...) do.call( p_LISA_method(name) , list(x@pointer , ... ))
})

# Constructors for p_LISA class
# Note: here simply using argc to determine which constructor should be called
#
setMethod( "initialize", "p_LISA", function(.Object, ...) {
  argv = list(...)
  argtypes <- mapply(class, argv);
  argc <- length(argtypes);

  if (argc == 0) {
    # this is for using p_LISA as a member in class('LISA')  in lisa.R
  } else {
    .Object@pointer <- argv[[1]]
  }
  .Object
})

