# rgeoda

rgeoda is a R package for spatial data analysis based on libgeoda and GeoDa. It provides spatial data analysis functionalities including Exploratory Spatial Data Analysis, Spatial Cluster Detection and Clustering Analysis, Regionalization, etc. based on the C++ source code of GeoDa, which is an open-source software tool that serves as an introduction to spatial data analysis. The GeoDa software and its documentation are available at https://geodacenter.github.io.
  
The rgeoda site is built using pkgdown: https://geodacenter.github.io/rgeoda


## Installation


```R
install.packages("rgeoda")
```

![cran status](https://www.r-pkg.org/badges/version/rgeoda)
![cran release](https://www.r-pkg.org/badges/last-release/rgeoda)
![cran downloads](https://cranlogs.r-pkg.org/badges/grand-total/rgeoda)

## Quick Start

```R
library(sf)
library(rgeoda)

guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry <- st_read(guerry_path)

w <- queen_weights(guerry)
lisa <- local_moran(w, guerry['Crm_prs'])
clusters <- skater(4, w, guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])
```
### Tutorials

https://geodacenter.github.io/rgeoda/articles/rgeoda_tutorial.html

### APIs

https://geodacenter.github.io/rgeoda/reference/
  
## Current version 0.0.8

* Map Classification
   * NaturalBreaks
   * QuantileBreaks
   * Hinge15Breaks
   * Hinge30Breaks
   * PercentileBreaks
   * StddevBreaks
   
* Spatial Weights
    * Queen
    * Rook
    * Distance based
    * K-Nearest Neighbor
    * Kernel
    * Read GAL/GWT/SWM weights
    
* Spatial Autocorrelation
    * Local Moran
    * Local Moran EB Rates
    * Local Geary
    * Local Getis-Ord 
    * Multivariate Local Geary
    * Local Join Count
    * Bivariate Local Join Count
    * (Multivariate) Colocation Local Join Count
    * Quantile LISA
    * Multivariate Quantile LISA
    * Neighbor Match Test

* Spatial Clustering
    * SCHC Spatial Constrained Hierarchical Clustering 
      * Single-linkage
      * Complete-linkage
      * Average-linkage
      * Ward-linkage
    * SKATER
    * REDCAP
      * First-order and Single-linkage
      * Full-order and Complete-linkage
      * Full-order and Average-linkage
      * Full-order and Single-linkage
      * Full-order and Ward-linkage
    * AZP
      * greedy
      * Tabu Search
      * Simulated Annealing
    * Max-p
      * greedy
      * Tabu Search
      * Simulated Annealing

## Build and install from source code

In R console, one can use devtools to install rgeoda from its **source package**:

```R
devtools::install_github("geodacenter/rgeoda")
```

#### Mac

For Mac users, the “Xcode Command Line Tools” need to be installed for installing rgeoda. It is a free software provided by Apple, which can be installed by using the following command in a terminal:
```
xcode-select --install 
```

Note that the Xcode tools are not automatically updated when a new version of Xcode is installed. In order to make
sure you have the latest version, use:

```
sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```

In order to make sure to have the correct C++ compiler for R 4.0 and later, follow the instructions
on https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/.


#### Windows

On Windows, the `Rtools` needs to be installed first. https://cran.r-project.org/bin/windows/Rtools/

#### Linux

For Linux users, the “Build Essential Tools” needs to be installed first.
```
sudo apt-get update
sudo apt-get install build-essential
```
