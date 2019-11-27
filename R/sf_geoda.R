random_string <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


geoda_to_sf = function(gda) {
  if (!require("sf")) {
    stop("package sf not available: install first?")
  }
}

#' @title Create a geoda object from a sf object
#' @description Create a geoda object from a sf object returned from 'st_read()' function
#' @param sf_obj  An instance of sf object
#' @param with_table  (Optional, Default: FALSE)If create a table from sf dataframe object.
#' @return geoda_obj An instance of geoda class
#' @export
sf_to_geoda = function(sf_obj, ...) {
  if (!require("sf")) {
    stop("package sf not available: install first?")
  }
  if (!require("wkb")) {
    stop("package wkb not available: install first?")
  }

  kwargs <- list(...)
  with_table <- ifelse(hasArg("with_table"), kwargs$with_table, TRUE)

  # geometries
  sf_geom <- st_geometry(sf_obj)
  geoms_wkb <- st_as_binary(sf_geom)
  n_obs <- length(sf_geom)
  wkb_bytes_len <- sapply(geoms_wkb, function(x){ return(length(x))})
  wkb_vec <- unlist(geoms_wkb)

  # in-memory name
  file_name <- random_string(1)

  # table
  sf_df <- as.data.frame(sf_obj)
  col_names <- colnames(sf_df)
  n_cols <- length(col_names)
  tbl <- GeoDaTable()
  if (with_table) {
    for (i in 1:n_cols) {
      col_nm <- col_names[[i]]
      if (col_nm == "geometry") next

      dat <- sf_df[, col_nm]
      ft <- class(dat)
      if (ft == "factor") {
        tbl$AddStringColumn(col_nm, dat)

      } else if (ft == "integer" || ft == "logical") {
        tbl$AddIntColumn(col_nm, dat)

      } else if (ft == "double" || ft == "numeric") {
        tbl$AddRealColumn(col_nm, dat)

      } else {
        dat <- as.character(dat)
        tbl$AddStringColumn(col_names[[i]], dat)
      }
    }

  } else {
    n_cols <- 0
    col_names <- rep("", 0)
  }
  # map_type
  map_type <- "map_polygons"
  geom_type <- st_geometry_type(sf_obj)[[1]]
  if (geom_type == "MULTIPOINT" || geom_type == "POINT") {
    map_type <- "map_points"
  } else if (geom_type == "MULTILINESTRING" || geom_type == "LINESTRING") {
    map_type <- "map_lines"
  }

  # prj4
  proj4_str <- st_crs(sf_obj)[[2]]
  gda <- GeoDa(file_name, map_type, n_obs, tbl, as.integer(wkb_vec), wkb_bytes_len, proj4_str)
  return(geoda$new(gda))
}

#' @title Create a geoda object from a sp object
#' @description The sp package has been an essential tool which provides spatial data-structures and many utility functions to do spatial analysis in R. It has been a core dependent library for many other packages, e.g. rgdal (IO), maptools (mapping), spdep (spatial weights, spatial statistics, and spatial models) etc.
#' Using rgdal to read a ESRI Shapefile will return a sp (Spatial object) object, which could be either a SpatialPointsDataFrame (using an AttributeList for its data slot directly), a SpatialLinesDataFrame, or a SpatialPolygonsDataFrame.
#' @param sp_obj  An instance of sp object
#' @param with_table  (Optional, Default: FALSE)If create a table from sp dataframe object.
#' @return geoda_obj An instance of GeoDa object
#' @export
sp_to_geoda = function(sp_obj, ...) {
  if (!require("sp")) {
    stop("package sp not available: install first?")
  }
  if (!require("wkb")) {
    stop("package wkb not available: install first?")
  }

  kwargs <- list(...)
  with_table <- ifelse(hasArg("with_table"), kwargs$with_table, FALSE)

  # geometries
  geoms_wkb <- writeWKB(sp_obj)
  n_obs <- length(geoms_wkb)
  wkb_bytes_len <- sapply(geoms_wkb, function(x) {return(length(x))})
  wkb_vec <- unlist(geoms_wkb)

  # in-memory name
  file_name <- random_string(1)

  # table
  col_names <- colnames(sp_obj@data)
  n_cols <- length(col_names)
  tbl <- GeoDaTable()
  if (with_table) {
    for (i in 1:n_cols) {
      ft <- class(sp_obj@data[[i]])
      if (ft == "factor") {
        dat <- sp_obj@data[[i]]
        tbl$AddStringColumn(col_names[[i]], dat)

      } else if (ft == "integer" || ft == "logical") {
        dat <- sp_obj@data[[i]]
        tbl$AddIntColumn(col_names[[i]], dat)

      } else if (ft == "double" || ft == "numeric") {
        dat <- sp_obj@data[[i]]
        tbl$AddRealColumn(col_names[[i]], dat)

      } else {
        dat <- as.character(sp_obj@data[[i]])
        tbl$AddStringColumn(col_names[[i]], dat)
      }
    }
  } else {
    n_cols <- 0
    col_names <- rep("", 0)
  }

  # map_type
  map_type <- "map_polygons"
  if (is(sp_obj, 'SpatialPointsDataFrame')) {
    map_type <- "map_points"
  } else if (is(sp_obj, "SpatialLinesDataFrame-class")) {
    map_type <- "map_lines"
  }
  # prj4
  proj4_str <- sp_obj@proj4string
  gda <- GeoDa(file_name, map_type, n_obs, tbl, as.integer(wkb_vec), wkb_bytes_len, proj4_str@projargs)
  return(geoda$new(gda))
}

#' create a sp object from a GeoDa object
#'
#' @export
geoda_to_sp = function(gda, geometry_only=TRUE) {
  if (!require("sp")) {
    stop("package sp not available: install first?")
  }
  if (!require("wkb")) {
    stop("package wkb not available: install first?")
  }
  # map_type
  map_type <- gda$GetMapType()

  # geometries
  n_obs <- gda$GetNumObs()
  wkb <- gda$GetGeometryWKB()
  wkb <- I(wkb)

  # create Spatial object
  sp_obj <- readWKB(
    wkb,
    #id = c("San Francisco", "New York"), using default sequential ids
    proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  )

  # create data frome
  df <- data.frame(row.names=row.names(sp_obj))
  if (geometry_only == FALSE) {
      n_cols <- gda$GetNumCols()
      col_nms <- gda$GetFieldNames()
      col_tps <- gda$GetFieldTypes()
      for (i in 1:n_cols) {
        col_nm <- col_nms[[i]]
        col_tp <- col_tps[[i]]
        if (col_tp == "integer") {
          df[, col_nm] <- gda$GetIntegerCol(col_nm)
        } else if (col_tp == "numeric") {
          df[, col_nm] <- gda$GetNumericCol(col_nm)
        } else {
          df[, col_nm] <- gda$GetStringCol(col_nm)
        }
      }
  }

  # create spatial dataframe
  if (map_type == "polygon_type") {
    return(SpatialPolygonsDataFrame(sp_obj, data = df))
  } else if(map_type == "point_type") {
    return(SpatialPointsDataFrame(sp_obj, data = df))
  } else if (map_type == "line_type") {
    return(SpatialLinesDataFrame(sp_obj, data = df))
  }
  return (NULL)
}
