#' @title A R wrapper for GeoDa
#' @description geoda is a RefClass that wraps the C++ GeoDa class (via _p_GeoDa defines in rgeoda.R)
#' @param o_gda An object of _p_GeoDa class
#' @return An object of geoda class
#' @export
geoda <- setRefClass("geoda",
  fields = list(
    gda = "_p_GeoDa",
    map_type = "character",
    n_cols = "numeric",
    n_obs = "numeric",
    field_names = "vector",
    field_types = "vector",
    table = "data.frame"
  ),
  methods = list(
    initialize = function(o_gda) {
      .self$gda = o_gda
      .self$map_type = gda$GetMapType()
      .self$n_cols = gda$GetNumCols()
      .self$n_obs = gda$GetNumObs()
      .self$field_names = gda$GetFieldNames()
      .self$field_types = gda$GetFieldTypes()
      .self$table <- data.frame()[1:.self$n_obs, ]

      if (length(.self$field_names) > 0) {
        for (i in 1:.self$n_cols) {
          f_nm <- .self$field_names[[i]]
          f_tp <- .self$field_types[[i]]
          if (f_tp == "numeric") {
            .self$table[[f_nm]] <- gda$GetNumericCol(f_nm)
          } else if (f_tp == "integer") {
            .self$table[[f_nm]] <- gda$GetIntegerCol(f_nm)
          } else {
            .self$table[[f_nm]] <- gda$GetStringCol(f_nm)
          }
        }
      }
    },
    GetNumCols = function(...) {
      return(gda$GetNumCols())
    },
    GetNumObs = function(...) {
      return(gda$GetNumObs())
    },
    GetFieldNames = function(...) {
      return(gda$GetFieldNames())
    },
    GetFieldTypes = function(...) {
      return(gda$GetFieldTypes())
    },
    GetMapType = function(...) {
      return(gda$GetMapType())
    },
    GetIntegerCol = function(col_name) {
      "Get the values (in integer type) from a column"
      return(gda$GetIntegerCol(col_name))
    },
    GetRealCol = function(col_name) {
      "Get the values (in real type) from a column"
      return(gda$GetNumericCol(col_name))
    },
    GetUndefinedVals = function(col_name) {
      return(gda$GetUndefinesCol(col_name))
    }
  )
)

#' @title Create a geoda object by reading a spatial dataset
#' @description Create a geoda object by reading a spatial dataset. The dataset that rgeoda supports includes: ESRI Shapefile, MapInfo File, CSV, GML, GPX, KML, GeoJSON, TopoJSON, OpenFileGDB, GFT Google Fusion Tables, CouchDB
#' @param ds_path (character) The path of the spatial dataset
#' @return gda_obj An object of geoda instance
#' @export
geoda_open <- function(ds_path) {
  if (typeof(ds_path) != "character") {
    stop ("Only a string of input datasource is allowed")
  }
  o_gda <- GeoDa(ds_path)
  return(geoda$new(o_gda))
}
