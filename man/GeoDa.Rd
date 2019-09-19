\name{GeoDa}
\alias{GeoDa}
\alias{sp_to_geoda}
\title{
  GeoDa()
}
\description{
There are two ways to construct a GeoDa object: the first one is passing the path of the spatial data as an argument; the second one is creating a GeoDa object from a sp/sf object. For latter, please refer to function sp_to_geoda() and sf_to_geoda().
}
\usage{
GeoDa(file_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path}{
  Path to the spatial data, e.g. ESRI Shapefile (.shp)
}
}
\details{

}
\value{
  The data formats that rgeoda can read directly:
  
  \item{}{ESRI Shapefile: ESRI Shapefile}
  \item{}{MapInfo File : MapInfo File}
  \item{}{CSV: Comma Separated Value (.csv)}
  \item{}{GML: Geography Markup Language (GML)}
  \item{}{GPX: GPX}
  \item{}{KML: Keyhole Markup Language (KML)}
  \item{}{GeoJSON: GeoJSON}
  \item{}{TopoJSON: TopoJSON}
  \item{}{OpenFileGDB: ESRI FileGDB}
  \item{}{GFT: Google Fusion Tables}
  \item{}{CouchDB: CouchDB / GeoCouch}
  \item{}{Carto: Carto}
}
\references{
https://lixun910.github.io/pygeoda/doc_intro.html
}
\author{
Xun Li (lixun910@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(rgeoda)

# In this example, we will use an ESRI Shapefile natregimes.shp comes with the package:
nat_path <- system.file("extdata", "natregimes.shp", package = "rgeoda")

# Use the path above (or the path to your own dataset), we can create a GeoDa instance, and the GeoDa instance is the main entry point of rgeoda.
gda <- GeoDa(nat_path)

# rgeoda provides some functions to check the meta data of the loaded dataset:
num_obs <- gda$GetNumObs()
num_cols <- gda$GetNumCols()
field_types <- gda$GetFieldTypes()
field_names <- gda$GetFieldNames()

# Access Table Data
# One can call function GetNumericCol(string col_name) to get numeric data (vector type) from GeoDa instance. For example, to get data of column “HR60”:
hr60 <- gda$GetNumericCol("HR60")

# Spatial Weights
# One can call function CreateContiguityWeights(string poly_id="", bool is_queen=true, int order=1, bool include_lower_order=false) to create a Queen or Rook contiguity weights:
queen_w <- gda$CreateContiguityWeights(is_queen=TRUE)

# To access the properties of the created weights object, one can just access the attributes of the returned weights object queen_w:

cat("weight_type: ", queen_w$weight_type, 
    "\nis_symmetric: ", queen_w$is_symmetric, 
    "\nsparsity:",queen_w$sparsity,
    "\ndensity:",queen_w$density, 
    "\nmin_nbrs:",queen_w$min_nbrs, 
    "\nmax_nbrs:",queen_w$max_nbrs,
    "\nmean_nbrs:",queen_w$mean_nbrs, 
    "\nmedian_nbrs:",queen_w$median_nbrs, "\n")
    
# Spatial Data Analysis
# Local Spatial Autocorrelation

# Using the created Queen weights queen_w and the data hr60, we can call function LISA(GeoDaWeight w, vector data) to compute the local spatial autocorrelation of variable “HR60”.
lisa <- gda$LISA(queen_w, hr60)

# We can access the LISA results by calling the “getter” methods from the returned LISA object:
lags <- lisa$GetLagValues()

# Get local moran values:
lms <- lisa$GetLocalMoranValues()

# Get pseudo-p values:
pvals <- lisa$GetLocalSignificanceValues()

# Get LISA category values:
cats <- lisa$GetClusterIndicators()

# You can easily re-run the LISA computation by calling its Run() function. For example, re-run the above LISA exampe using 9999 permutations:
lisa$SetNumPermutations(9999)
lisa$Run()

# Since rgeoda is using GeoDa’s C++ code, by default, rgeoda uses multi-threading to accelerate the computation of Local Moran. One can also specify how many threads to run the computation:
lisa$SetNumThreads(4)
lisa$Run()

# Clustering

# SKATER
# Using the function SKATER(int k, GeoDaWeight w, Vector column_names, String distance_method='euclidean') to run a spatially constrained clustering on current dataset.
skater <- gda$SKATER(10, queen_w, c("HR60", "PO60"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GeoDa }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~LISA }% __ONLY ONE__ keyword per line
