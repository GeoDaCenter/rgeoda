############################################################
#' @title Spatial C(K)luster Analysis by Tree Edge Removal
#' @description SKATER forms clusters by spatially partitioning data that has similar values for features of interest.
#' @param k The number of clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param bound_vals (optional) A 1-d vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (int,optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- skater(4, queen_w, data)
#' guerry_clusters
#' @export
skater <- function(k, w, data, bound_vals=vector('numeric'), min_bound=0, distance_method="euclidean", random_seed=123456789, cpu_threads=6) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is less than the number of observations.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }
  return(p_skater(k, w$GetPointer(), data, distance_method, bound_vals, min_bound, random_seed, cpu_threads))
}


############################################################
#' @title Spatially Constrained Hierarchical Clucstering (SCHC)
#' @description Spatially constrained hierarchical clustering is a special form of constrained clustering, where the constraint is based on contiguity (common borders).
#' The method builds up the clusters using agglomerative hierarchical clustering methods:
#' single linkage, complete linkage, average linkage and Ward's method (a special form of centroid linkage).
#' Meanwhile, it also maintains the spatial contiguity when merging two clusters.
#' @param k The number of clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param method {"single", "complete", "average","ward"}
#' @param bound_vals (optional) A 1-d vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- schc(4, queen_w, data, "complete")
#' guerry_clusters
#' @export
schc <- function(k, w, data, method="average", bound_vals=vector('numeric'), min_bound=0, distance_method="euclidean") {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is less than the number of observations.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  method_cands <- c("single", "complete", "average","ward")
  if (!(method %in% method_cands)) {
    stop("The SCHC method has to be one of ", method_cands)
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }
  return(p_schc(k, w$GetPointer(), data, method, distance_method, bound_vals, min_bound))
}


############################################################
#' @title Regionalization with dynamically constrained agglomerative clustering and partitioning
#' @description REDCAP (Regionalization with dynamically constrained agglomerative
#' clustering and partitioning) is developed by D. Guo (2008). Like SKATER, REDCAP
#' starts from building a spanning tree with 4 different ways (single-linkage, average-linkage, ward-linkage
#' and the complete-linkage). The single-linkage way leads to build a minimum spanning tree.
#' Then,REDCAP provides 2 different ways (first-order and full-order constraining) to
#' prune the tree to find clusters. The first-order approach with a minimum spanning tree is
#' exactly the same with SKATER. In GeoDa and pygeoda, the following methods are provided:
#' \* First-order and Single-linkage
#' \* Full-order and Complete-linkage
#' \* Full-order and Average-linkage
#' \* Full-order and Single-linkage
#' \* Full-order and Ward-linkage
#' @param k The number of clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param method {"firstorder-singlelinkage", "fullorder-completelinkage", "fullorder-averagelinkage","fullorder-singlelinkage", "fullorder-wardlinkage"}
#' @param bound_vals (optional) A 1-d vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (int,optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- redcap(4, queen_w, data, "fullorder-completelinkage")
#' guerry_clusters
#' @export
redcap <- function(k, w, data, method="fullorder-averagelinkage", bound_vals=vector('numeric'), min_bound=0, distance_method="euclidean", random_seed=123456789, cpu_threads=6) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is less than the number of observations.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  method_cands <- c("firstorder-singlelinkage", "fullorder-completelinkage", "fullorder-averagelinkage","fullorder-singlelinkage", "fullorder-wardlinkage")
  if (!(method %in% method_cands)) {
    stop("The REDCAP method has to be one of ", method_cands)
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }
  return(p_redcap(k, w$GetPointer(), data, method, distance_method, bound_vals, min_bound, random_seed, cpu_threads))
}

############################################################
#' @title A greedy algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained clustering where a finite number of geographical areas, n, are aggregated into the maximum number of regions, p, such that each region satisfies the following const raints: 1. The areas within a region must be geographically connected.
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param bound_vals A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param iterations (optional): The number of iterations of greedy algorithm. Defaults to 99.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_vals <- guerry_df['Pop1831'][,1]
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_greedy(queen_w, data, bound_vals, min_bound, iterations=99)
#' maxp_clusters
#' }
#' @export
maxp_greedy <- function(w, data, bound_vals, min_bound, iterations=99, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789, cpu_threads=6) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (length(bound_vals) != w$num_obs) {
    stop("The bound_vals has to be a list of numeric values, e.g. a column of input table.")
  }
  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_maxp_greedy(w$GetPointer(), data, bound_vals, min_bound, iterations, initial_regions, distance_method, random_seed, cpu_threads))
}

############################################################
#' @title A simulated annealing algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained clustering where a finite number of geographical areas, n, are aggregated into the maximum number of regions, p, such that each region satisfies the following const raints: 1. The areas within a region must be geographically connected.
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param bound_vals A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param cooling_rate The cooling rate of a simulated annealing algorithm. Defaults to 0.85
#' @param iterations (optional): The number of iterations of SA algorithm. Defaults to 99.
#' @param sa_maxit (optional): The number of iterations of simulated annealing. Defaults to 1
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_vals <- guerry_df['Pop1831'][,1]
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_sa(queen_w, data, bound_vals, min_bound, cooling_rate=0.85, sa_maxit=1)
#' maxp_clusters
#' }
#' @export
maxp_sa <- function(w, data, bound_vals, min_bound, cooling_rate, sa_maxit=1, iterations=99, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789, cpu_threads=6) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (length(bound_vals) != w$num_obs) {
    stop("The bound_vals has to be a list of numeric values, e.g. a column of input table.")
  }
  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_maxp_sa(w$GetPointer(), data, bound_vals, min_bound, iterations, cooling_rate, sa_maxit, initial_regions, distance_method, random_seed, cpu_threads))
}

############################################################
#' @title A tabu-search algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained clustering where a finite number of geographical areas, n, are aggregated into the maximum number of regions, p, such that each region satisfies the following const raints: 1. The areas within a region must be geographically connected.
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param bound_vals A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param tabu_length (optional): The length of a tabu search heuristic of tabu algorithm. Defaults to 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param iterations (optional): The number of iterations of Tabu algorithm. Defaults to 99.
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_vals <- guerry_df['Pop1831'][,1]
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_tabu(queen_w, data, bound_vals, min_bound, tabu_length=10, conv_tabu=10)
#' maxp_clusters
#' }
#' @export
maxp_tabu <- function(w, data, bound_vals, min_bound, tabu_length=10, conv_tabu=10, iterations=99, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789, cpu_threads=6) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (length(bound_vals) != w$num_obs) {
    stop("The bound_vals has to be a list of numeric values, e.g. a column of input table.")
  }
  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_maxp_tabu(w$GetPointer(), data, bound_vals, min_bound, iterations, tabu_length, conv_tabu, initial_regions, distance_method, random_seed, cpu_threads))
}

############################################################
#' @title A greedy algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand, as in most other clustering methods considered so far.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param bound_vals (optional) A numeric vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_greedy(5, queen_w, data)
#' azp_clusters
#' }
#' @export
azp_greedy <- function(p, w, data, bound_vals=vector('numeric'), min_bound=0, inits=0, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_azp_greedy(p, w$GetPointer(), data, bound_vals, min_bound, inits, initial_regions, distance_method, random_seed))
}

############################################################
#' @title A simulated annealing algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand, as in most other clustering methods considered so far.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param cooling_rate The cooling rate of a simulated annealing algorithm. Defaults to 0.85
#' @param sa_maxit (optional): The number of iterations of simulated annealing. Defaults to 1
#' @param bound_vals (optional) A numeric vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_sa(5, queen_w, data, cooling_rate = 0.85)
#' azp_clusters
#' }
#' @export
azp_sa<- function(p, w, data, cooling_rate, sa_maxit=1, bound_vals=vector('numeric'), min_bound=0, inits=0, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_azp_sa(p, w$GetPointer(), data, cooling_rate, sa_maxit, bound_vals, min_bound, inits, initial_regions, distance_method, random_seed))
}

############################################################
#' @title A tabu algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand, as in most other clustering methods considered so far.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param data A list of numeric vectors of selected variable
#' @param tabu_length The length of a tabu search heuristic of tabu algorithm. e.g. 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param bound_vals (optional) A numeric vector of selected bounding variable
#' @param min_bound (optional) A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @return A list of numeric vectors represents a group of clusters
#' @examples
#' \dontrun{
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- geoda_open(guerry_path)
#' queen_w <- queen_weights(guerry)
#' guerry_df <- as.data.frame(guerry) # use as data.frame
#' data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_tabu(5, queen_w, data, tabu_length=10, conv_tabu=10)
#' azp_clusters
#' }
#' @export
azp_tabu<- function(p, w, data, tabu_length=10, conv_tabu=10, bound_vals=vector('numeric'), min_bound=0, inits=0, initial_regions=vector('numeric'), distance_method="euclidean", random_seed=123456789) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (length(data) < 1) {
    stop("The data from selected variable is empty.")
  }
  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  return(p_azp_tabu(p, w$GetPointer(), data, tabu_length, conv_tabu, bound_vals, min_bound, inits, initial_regions, distance_method, random_seed))
}


############################################################
#' @title Between Sum of Square
#' @description Compute between sum of square value of a group of clusters
#' @param clusters A list of numeric vectors which returns from spatial clustering methods, e.g. skater()
#' @param data A list of numeric vectors which is used in spatial clustering methods, e.g. skater(k, w, data)
#' @return A value of between sum of square value
#' @export
between_sumofsquare <- function(clusters, data) {
  return(p_betweensumofsquare(clusters, data))
}


############################################################
#' @title Total Sum of Square
#' @description Compute total sum of square value of a 2d tuple data
#' @param data A list of numeric vectors which is used in spatial clustering methods, e.g. skater(k, w, data)
#' @return A value of total sum of square value
#' @export
total_sumofsquare <- function(data) {
  return(p_totalsumofsquare(data))
}

############################################################
#' @title Within Sum of Square
#' @description Compute within sum of square value of a group of clusters
#' @param clusters A list of numeric vectors which returns from spatial clustering methods, e.g. skater()
#' @param data A list of numeric vectors which is used in spatial clustering methods, e.g. skater(k, w, data)
#' @return A value of within sum of square value
#' @export
within_sumofsquare <- function(clusters, data) {
  return(p_withinsumofsquare(clusters, data))
}

