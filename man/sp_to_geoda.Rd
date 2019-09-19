\name{sp_to_geoda}
\alias{sp_to_geoda}
\alias{GeoDa}
\title{
  sp_to_geoda()
}
\description{
This function creates a rgeoda object from the sp object.

The sp package has been an essential tool which provides spatial data-structures and many utility functions to do spatial analysis in R. It has been a core dependent library for many other packages, e.g. rgdal (IO), maptools (mapping), spdep (spatial weights, spatial statistics, and spatial models) etc.

Using rgdal to read a ESRI Shapefile will return a sp (Spatial object) object, which could be either a SpatialPointsDataFrame (using an AttributeList for its data slot directly), a SpatialLinesDataFrame, or a SpatialPolygonsDataFrame.

}
\usage{
sp_to_geoda(sp_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp_obj}{
An instance of sp object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
  \item{geoda}{An instance of GeoDa}
%% ...
}
\references{
https://lixun910.github.io/pygeoda/wn1_rgeoda_sp.html
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
guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")

# Create a sp object from the Guerry.shp file
library(rgdal)
guerry_sp <- readOGR(guerry_path)

# Create a rgeoda object from the sp object guerry_sp you just created.
guerry_gda <- sp_to_geoda(guerry_sp)

# That’s it! It should be very fast: rgeoda uses wkb, which is a binary representation of geometries, to exchange data between sp and rgeoda in memory.

# Now, with the rgeoda object guerry_gda, you can call GeoDa’s spatial analysis functions easily. For example, we want to examine the local Moran of variable “crm_prs” (Population per Crime against persons):
queen_w <- guerry_gda$CreateContiguityWeights()
crm_lisa <- guerry_gda$LISA(queen_w, as.numeric(paste(guerry_sp$Crm_prs)))
crm_lisa_clusters <- crm_lisa$GetClusterIndicators()

# Now, with the LISA results, we can do exploratory spatial data analysis by generating a LISA cluster map:
plot(guerry_sp, 
     col=sapply(crm_lisa_clusters, function(x){return(lisa_colors[[x+1]])}), 
     border = "#333333", lwd=0.2)
title(main = "Local Moran Map")
legend('bottomright', legend = lisa_labels, fill = lisa_colors, border = "#eeeeee")

# From the above code, you can see that we still use sp object to do plotting. The values of cluster indicators from rgeoda’s LISA object are used to make the LISA map.This emphasis that rgeoda is a API focused library, which only provides the core functions of spatial data analysis for easy integration.

# You can easily append the lisa results to original sp object by manipulating the data.frame object inside the sp object.

# If you check the values of the cluster indicators, you will see they are integer numbers 0 (not significant), 1 (high-high cluster), 2 (low-low cluster), 3 (low-high cluster), 4 (high-low cluster), 5 (neighborless/island), 6 (undefined), which are excatly the same with GeoDa software when you save LISA results to a table:
crm_lisa_clusters

# We also specify what color is associated with which cluster value to replace the lisa_colors and lisa_labels, which are defined in sf_geoda.R:
lisa_colors <- c("#eeeeee","#FF0000","#0000FF","#a7adf9", "#f4ada8", "#999999")

# We then create an array of hex color, each observation is set to a hex color according to its cluster indicator value.
sapply(crm_lisa_clusters, function(x){return(lisa_colors[[x+1]])})


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sp }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~ESDA }% __ONLY ONE__ keyword per line

