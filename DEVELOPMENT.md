## For developers:

### Prerequisite
pre-requisite for dev work on OSX 10.13:

1. Installing Xcode Command Line Tools 
xcode-select --install 

2. Install SWIG
brew install swig

### Add new features for rGeoDa:

If you are interested in adding new features in rgeoda, here is what you can do to build/run (later write your code and debug) libgeoda-r from source:

1. Use another repo `libgeoda`,  git clone https://github.com/lixun910/libgeoda
2. cd libgeoda/swig/R

To build rgeoda for development:
```
./test.sh
```
You should get three files generated:
```
rgeoda.R
rgeoda.cpp
rgeoda.so
```

With the three files, you can start dev work. Check the file test.R
```
dyn.load("rgeoda.so")
source("rgeoda.R")
source("sf_geoda.R")
cacheMetaData(1)

gda <- GeoDa("../../data/Guerry.shp")
```

The functions that I am working on are in the file `sf_geoda.R`. I added three functions, `sf_to_geoda`,  `sp_to_geoda` and `geoda_to_sp`. The plan is to add one more function `geoda_to_sf` to allow rgeoda working with sf (simple features) seamlessly.

### Logic in `sf_geoda.R`

If you read the `sf_geoda.R`, you will see how to call functions from rgeoda. But here are some useful functions that you might want to know when work with sf:

TWO scenarios using rgeoda:  

1. user load spatial data using sf/sp, then create a`GeoDa` instance to call functions of spatial analysis 
2. user load spatial data to create a `GeoDa` instance, then create a sf/sp instance from `GeoDa` instance.

So, there are two ways to create an instance of `GeoDa` in R
```R
GeoDa("path/to/spatial_data_file")
```
>  ESRI Shapefile -vector- (rw+v): ESRI Shapefile
 MapInfo File -vector- (rw+v): MapInfo File
 CSV -vector- (rw+v): Comma Separated Value (.csv)
 GML -vector- (rw+v): Geography Markup Language (GML)
 GPX -vector- (rw+v): GPX
 KML -vector- (rw+v): Keyhole Markup Language (KML)
 GeoJSON -vector- (rw+v): GeoJSON
 TopoJSON -vector- (rov): TopoJSON
 OpenFileGDB -vector- (rov): ESRI FileGDB
 GFT -vector- (rw+): Google Fusion Tables
 CouchDB -vector- (rw+): CouchDB / GeoCouch
 Carto -vector- (rw+): Carto

```R
GeoDa(
    layer_name : character,  
    map_type  : character, 
    num_obs : integer,  
    geoda_table : GeoDaTable, 
    wkb : integer,
    wkb_bytes_length : array (integer)
    proj4_str : character
)
```

####  Scenario 1

Create a `sp` object (a spatial datframe) from a `GeoDa` object.

The spatial dataframe could be one of the 3 types: `SpatialPolygonsDataFrame`, `SpatialPointsDataFrame`, `SpatialLinesDataFrame`. 

There are (at least) two input parameters needed to create an object of spatial dataframe: a `spatial object` and a `data.frame` object.

To create a `spatial object`, we call `GeoDa$GetGeometryWKB()` to get the geometries (as a list of wkb objects) from a `geoda` object.  
```R
wkb <- gda$GetGeometryWKB()
wkb <- I(wkb)
```
> I() function is to isolates or insulates the contents of wkb list, that is required by `readWKB()` function

Then, you can use the list of wkb object to create a spatial object:
```R
  sp_obj <- readWKB(
    wkb,
    #id = c("San Francisco", "New York"), using default sequential ids
    proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  )
```

I think this will be the same logic to create a `sf` object using `wkb`.

To create a R data.frame,  you need to loop the fields/columns in the `GeoDa` object. You can ref the code block in `geoda_to_sp` function:
```R
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
```

The final step is to create a spatial dataframe ( I think sf object should have different forms, but still be an extension of data.frame):

```R
# create spatial dataframe
  if (map_type == "polygon_type") {
    return(SpatialPolygonsDataFrame(sp_obj, data = df))
  } else if(map_type == "point_type") {
    return(SpatialPointsDataFrame(sp_obj, data = df))
  } else if (map_type == "line_type") {
    return(SpatialLinesDataFrame(sp_obj, data = df))
  }
```

#### Scenario 2

In this scenario, `libgeoda` will create a in-memory dataset using the content of a `sp` object (the same design is applied to `sf` object). So, there is a random string for "file_name" to mark this in-memory dataset. Here are the details of other parameters:

##### Geometries

`wkb` package is used to create a raw wkb object (to represent the geometries) from a `sp` object. 
```R
 geoms_wkb <- writeWKB(sp_obj)
```
You will see `geoms_wkb` is a list of raw bytes array. To allow libgeoda accessing its content, we need to flatten this list
```R
wkb_vec <- unlist(geoms_wkb)
```
Besides, we also need to tell libgeoda the size of the list, and the size of each raw array.
```R
n_obs <- length(geoms_wkb)
wkb_bytes_len <- sapply(geoms_wkb, function(x) {return(length(x))})
```
##### Attributes/Table

For attributes data, libgeoda will create an instance of `GeoDaTable()` with three member functions:
```
AddStringColumn(column_name : character, data : array of factor/chars)

AddIntColumn(column_name : character, data : array of integer)

AddRealColumn(column_name : character, data : array of numeric)
```
In sp_to_geoda function, you will see the code block that loop the columns in data.frame of a `sp` object, and create a `GeoDaTable` instance (the same logic should apply to sf_to_geoda):
```R
col_names <- colnames(sp_obj@data)
  n_cols <- length(col_names)
  tbl <- GeoDaTable()
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
```
##### MapType

`map_type` is hard-code character, and the values are:
map_polygons (ref to SpatialPolygonsDataFrame)
map_points (ref to SpatialPointsDataFrame)
map_lines (ref to SpatialLinesDataFrame)

##### Projection

`proj4_str`: projection information coded in proj4 format
```R
proj4_str <- sp_obj@proj4string
```


## APIs:

```R
GeoDa("path/to/spatial_data_file")
```

```R
CreateContiguityWeights(
           is_queen : TRUE or FALSE (TRUE default)
           polyid : character (default: "")
           order : integer (default =1
          include_lower_order : Logical (default: FALSE)
          precision_threshold : numeric (default: 0)
```

```R
LISA(GeoDaWeight, data : numeric list)
```

