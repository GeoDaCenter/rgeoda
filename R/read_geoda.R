#' @title A R wrapper for p_GeoDa
#' @description geoda is a RefClass that wraps the C++ GeoDa class (via p_GeoDa defines in rgeoda.R)
#' @field gda A pointer to C++ GeoDa object
#' @field map_type The map type: Point, Polygon (or LineSegment)
#' @field n_cols The number of columns
#' @field n_obs The number of observations
#' @field field_names A list of field names
#' @field field_types A list of field types  (integer, real, string)
#' @export
geoda <- setRefClass("geoda",
  fields = list(
    gda = "p_GeoDa",
    map_type = "integer",
    n_cols = "numeric",
    n_obs = "numeric",
    field_names = "vector",
    field_types = "vector",
    table = "data.frame"
  ),
  methods = list(
    initialize = function(o_gda) {
      "Constructor with a geoda object (internally used)"
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
          if (f_tp == "real") {
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
      "Get the number of columns"
      return(gda$GetNumCols())
    },
    GetNumObs = function(...) {
      "Get the number of observations"
      return(gda$GetNumObs())
    },
    GetFieldNames = function(...) {
      "Get the field names of all columns"
      return(gda$GetFieldNames())
    },
    GetFieldTypes = function(...) {
      "Get the field types (integer, real, string) of all columns"
      return(gda$GetFieldTypes())
    },
    GetMapType = function(...) {
      "Get the map type"
      return(gda$GetMapType())
    },
    GetIntegerCol = function(col_name) {
      "Get the integer values from a column"
      return(gda$GetIntegerCol(col_name))
    },
    GetRealCol = function(col_name) {
      "Get the real values from a column"
      return(gda$GetNumericCol(col_name))
    },
    GetUndefinedVals = function(col_name) {
      "Get the undefined flags from a column"
      return(gda$GetNullValues(col_name))
    },
    GetPointer = function() {
      "Get the C++ object pointer (internally used)"
      return(gda$GetPointer())
    }
  )
)

#' @title Create a geoda object by reading a spatial dataset
#' @description Create a geoda object by reading a spatial dataset. The dataset that rgeoda supports includes: ESRI Shapefile, MapInfo File, CSV, GML, GPX, KML, GeoJSON, TopoJSON, OpenFileGDB, GFT Google Fusion Tables, CouchDB
#' @param ds_path (character) The path of the spatial dataset
#' @return gda_obj An object of geoda instance
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' guerry_df <- as.data.frame(guerry) # access as a data.frame
#' head(guerry_df)
#' @export
geoda_open <- function(ds_path) {
  if (typeof(ds_path) != "character") {
    stop ("Only a string of input datasource is allowed")
  }
  o_gda <- p_GeoDa(ds_path)
  return(geoda$new(o_gda))
}

#' @title convert rgeoda instance to data.frame
#' @description Override the as.data.frame function for rgeoda instance
#' @param x A rgeoda object
#' @param row.names NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional optional parameters
#' @param ... other arguments passed to methods
#' @return A data.frame object
#' @export
as.data.frame.geoda <- function(x, row.names = NULL, optional = FALSE, ...) {
  return (x$table)
}
