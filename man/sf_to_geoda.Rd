\name{sf_to_geoda}
\alias{sf_to_geoda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
sf_to_geoda()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function creates a rgeoda object from the sf object.

The sf package has been popular tool to handle geospatial data. It is a good substitue of sp package which will be deprecated soon. sf package makes it easy to work with geospatial data. It also provides functions to do basic spatial data analysis. rgeoda provides helper function sf_to_geoda to create a GeoDa instance from a sf object.

}
\usage{
sf_to_geoda(sf_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf_obj}{
An instance of sf object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{geoda}{An instance of GeoDa}
}
\references{
https://lixun910.github.io/rgeoda/esda_rgeoda_sf.html
}
\author{
%%  ~~who you are~~
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

# Create a sf object from the Guerry.shp file
library(sf)
guerry_sf <- st_read(guerry_path)

# You can simply call plot() function to render the first 9 chorepleth maps using the frist 9 variables in the dataset:
plot(guerry_sf)

# rgeoda provides helper function sf_to_geoda to create a GeoDa instance from a sf object. rgeoda uses wkb, which is a binary representation of geometries, to exchange data between sf and libgeoda in memory.
guerry_gda <- sf_to_geoda(guerry_sf)

# with the rgeoda object guerry_gda, you can call GeoDa’s spatial analysis functions. For example, to examine the local Moran of variable “crm_prs” (Population per Crime against persons):
queen_w <- guerry_gda$CreateContiguityWeights()
crm_lisa <- guerry_gda$LISA(queen_w, as.numeric(paste(guerry_sf$Crm_prs)))
crm_lisa_clusters <- crm_lisa$GetClusterIndicators()

# Now, with the LISA results, we can do exploratory spatial data analysis by generating a LISA cluster map:
colors <- sapply(crm_lisa_clusters, function(x){return(lisa_colors[[x+1]])})
plot(st_geometry(guerry_sf), col = colors, border = "#333333", lwd=0.2)
title(main = "LISA of crm_prs")
legend('bottomright', legend = lisa_labels, 
       fill = lisa_colors, border = "#eeeeee")

# From the above code, you can see that we still use sf object to do plotting. The values of cluster indicators from rgeoda’s LISA object are used to make the LISA map.
# You can easily append the lisa results to original sf object by manipulating the data.frame object inside the sf object.
# We then create an array of hex color, each observation is set to a hex color according to its cluster indicator value.
sapply(crm_lisa_clusters, function(x){return(lisa_colors[[x+1]])})

# If you check the values of the cluster indicators, you will see they are integer numbers 0 (not significant), 1 (high-high cluster), 2 (low-low cluster), 3 (low-high cluster), 4 (high-low cluster), 5 (neighborless/island), 6 (undefined), which are excatly the same with GeoDa software when you save LISA results to a table:
crm_lisa_clusters

# We also specify what color is associated with which cluster value to replace the lisa_colors and lisa_labels, which are defined in sf_geoda.R:
lisa_colors <- c("#eeeeee","#FF0000","#0000FF","#a7adf9", "#f4ada8", "#999999")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sf }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~ESDA }% __ONLY ONE__ keyword per line
