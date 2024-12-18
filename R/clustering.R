############################################################
#' @title Spatial C(K)luster Analysis by Tree Edge Removal
#' @description SKATER forms clusters by spatially partitioning data that has
#' similar values for features of interest.
#' @param k The number of clusters
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only.
#' E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param bound_variable (optional) A data frame with selected bound variable
#' @param min_bound (optional) A minimum bound value that applies to all
#' clusters
#' @param scale_method One of the scaling methods ('raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the
#' distance betwen observation i and j. Defaults to "euclidean". Options are
#' "euclidean" and "manhattan"
#' @param random_seed (int,optional) The seed for random number generator.
#' Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel
#' computation
#' @param rdist (optional) The distance matrix (lower triangular matrix,
#' column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares",
#' "Within-cluster sum of squares", "Total within-cluster sum of squares",
#' and "The ratio of between to total sum of squares".
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- skater(4, queen_w, data)
#' guerry_clusters
#' @export
skater <- function(k, w, df, bound_variable=data.frame(), min_bound=0,
                   scale_method="standardize", distance_method="euclidean",
                   random_seed=123456789, cpu_threads=6, rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is
         less than the number of observations.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize',
                     'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean',
         'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }
  return(p_skater(k, w$GetPointer(), df, n_vars, scale_method, distance_method,
                  bound_values, min_bound, random_seed, cpu_threads, rdist))
}


############################################################
#' @title Spatially Constrained Hierarchical Clucstering (SCHC)
#' @description Spatially constrained hierarchical clustering is a special form of constrained clustering, where the constraint is based on contiguity (common borders).
#' The method builds up the clusters using agglomerative hierarchical clustering methods:
#' single linkage, complete linkage, average linkage and Ward's method (a special form of centroid linkage).
#' Meanwhile, it also maintains the spatial contiguity when merging two clusters.
#' @param k The number of clusters
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param method {"single", "complete", "average","ward"}
#' @param bound_variable (optional) A data frame with selected bound variabl
#' @param min_bound (optional) A minimum bound value that applies to all clusters
#' @param scale_method One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- schc(4, queen_w, data, "complete")
#' guerry_clusters
#' @export
schc <- function(k, w, df, method="average", bound_variable=data.frame(), min_bound=0, scale_method="standardize", distance_method="euclidean", rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is less than the number of observations.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')")
  }

  method_cands <- c("single", "complete", "average","ward")
  if (!(method %in% method_cands)) {
    stop("The SCHC method has to be one of {'single', 'complete', 'average','ward'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_schc(k, w$GetPointer(), df, n_vars, scale_method, method, distance_method, bound_values, min_bound, rdist))
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
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param method {"firstorder-singlelinkage", "fullorder-completelinkage", "fullorder-averagelinkage","fullorder-singlelinkage", "fullorder-wardlinkage"}
#' @param bound_variable (optional) A data frame with selected bound variabl
#' @param min_bound (optional) A minimum bound value that applies to all clusters
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (int,optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' guerry_clusters <- redcap(4, queen_w, data, "fullorder-completelinkage")
#' guerry_clusters
#' }
#' @export
redcap <- function(k, w, df, method="fullorder-averagelinkage", bound_variable=data.frame(), min_bound=0, scale_method="standardize", distance_method="euclidean", random_seed=123456789, cpu_threads=6, rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (k <1 && k > w$num_obs) {
    stop("The number of clusters should be a positive integer number, which is less than the number of observations.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  method_cands <- c("firstorder-singlelinkage", "fullorder-completelinkage", "fullorder-averagelinkage","fullorder-singlelinkage", "fullorder-wardlinkage")
  if (!(method %in% method_cands)) {
    stop("The REDCAP method has to be one of {'firstorder-singlelinkage', 'fullorder-completelinkage', 'fullorder-averagelinkage','fullorder-singlelinkage', 'fullorder-wardlinkage'}")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_redcap(k, w$GetPointer(), df, n_vars, method, scale_method, distance_method, bound_values, min_bound, random_seed, cpu_threads, rdist))
}

############################################################
#' @title A greedy algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained
#' clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is
#' geographically connected and the clusters could maximize internal homogeneity.
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param bound_variable A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param iterations (optional): The number of iterations of greedy algorithm. Defaults to 99.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_variable <- guerry['Pop1831']
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_greedy(queen_w, data, bound_variable, min_bound, iterations=99)
#' maxp_clusters
#' }
#' @export
maxp_greedy <- function(w, df, bound_variable, min_bound, iterations=99, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, cpu_threads=6, rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  if (inherits(bound_variable, "data.frame") == FALSE) {
    stop("The bound_variable needs to be a data.frame.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_maxp_greedy(w$GetPointer(), df, n_vars, bound_values, min_bound, iterations, initial_regions, scale_method, distance_method, random_seed, cpu_threads, rdist))
}

############################################################
#' @title A simulated annealing algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained
#' clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is
#' geographically connected and the clusters could maximize internal homogeneity.
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param bound_variable A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param cooling_rate The cooling rate of a simulated annealing algorithm. Defaults to 0.85
#' @param iterations (optional): The number of iterations of SA algorithm. Defaults to 99.
#' @param sa_maxit (optional): The number of iterations of simulated annealing. Defaults to 1
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_variable <- guerry['Pop1831']
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_sa(queen_w, data, bound_variable, min_bound, cooling_rate=0.85, sa_maxit=1)
#' maxp_clusters
#' }
#' @export
maxp_sa <- function(w, df, bound_variable, min_bound, cooling_rate, sa_maxit=1, iterations=99, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, cpu_threads=6, rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  if (inherits(bound_variable, "data.frame") == FALSE) {
    stop("The bound_variable needs to be a data.frame.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_maxp_sa(w$GetPointer(), df, n_vars, bound_values, min_bound, iterations, cooling_rate, sa_maxit, initial_regions, scale_method, distance_method, random_seed, cpu_threads, rdist))
}

############################################################
#' @title A tabu-search algorithm to solve the max-p-region problem
#' @description The max-p-region problem is a special case of constrained
#' clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is
#' geographically connected and the clusters could maximize internal homogeneity.
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param bound_variable A numeric vector of selected bounding variable
#' @param min_bound A minimum value that the sum value of bounding variable int each cluster should be greater than
#' @param tabu_length (optional): The length of a tabu search heuristic of tabu algorithm. Defaults to 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param iterations (optional): The number of iterations of Tabu algorithm. Defaults to 99.
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param initial_regions (optional): The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' bound_variable <- guerry['Pop1831']
#' min_bound <- 3236.67 # 10% of Pop1831
#' maxp_clusters <- maxp_tabu(queen_w, data, bound_variable, min_bound, tabu_length=10, conv_tabu=10)
#' maxp_clusters
#' }
#' @export
maxp_tabu <- function(w, df, bound_variable, min_bound, tabu_length=10, conv_tabu=10, iterations=99, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, cpu_threads=6, rdist=numeric()) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  if (min_bound <= 0) {
    stop("The min_bound has to be a positive numeric value.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  if (inherits(bound_variable, "data.frame") == FALSE) {
    stop("The bound_variable needs to be a data.frame.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_maxp_tabu(w$GetPointer(), df, n_vars, bound_values, min_bound, iterations, tabu_length, conv_tabu, initial_regions, scale_method, distance_method, random_seed, cpu_threads, rdist))
}

############################################################
#' @title A greedy algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param bound_variable (optional) A data frame with selected bound variabl
#' @param min_bound (optional) A minimum bound value that applies to all clusters
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_greedy(5, queen_w, data)
#' azp_clusters
#' }
#' @export
azp_greedy <- function(p, w, df, bound_variable=data.frame(), min_bound=0, inits=0, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, rdist=numeric()) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_azp_greedy(p, w$GetPointer(), df, n_vars, bound_values, min_bound, inits, initial_regions, scale_method, distance_method, random_seed, rdist))
}

############################################################
#' @title A simulated annealing algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param cooling_rate The cooling rate of a simulated annealing algorithm. Defaults to 0.85
#' @param sa_maxit (optional): The number of iterations of simulated annealing. Defaults to 1
#' @param bound_variable (optional) A data frame with selected bound variabl
#' @param min_bound (optional) A minimum bound value that applies to all clusters
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_sa(5, queen_w, data, cooling_rate = 0.85)
#' azp_clusters
#' }
#' @export
azp_sa<- function(p, w, df, cooling_rate, sa_maxit=1, bound_variable=data.frame(), min_bound=0, inits=0, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, rdist=numeric()) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_azp_sa(p, w$GetPointer(), df, n_vars, cooling_rate, sa_maxit, bound_values, min_bound, inits, initial_regions, scale_method, distance_method, random_seed, rdist))
}

############################################################
#' @title A tabu algorithm to solve the AZP problem
#' @description The automatic zoning procedure (AZP) was initially outlined in Openshaw (1977) as a way to address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares as a criterion of homogeneity. The number of regions needs to be specified beforehand.
#' @param p The number of spatially constrained clusters
#' @param w An instance of Weight class
#' @param df A data frame with selected variables only. E.g. guerry[c("Crm_prs", "Crm_prp", "Litercy")]
#' @param tabu_length The length of a tabu search heuristic of tabu algorithm. e.g. 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param bound_variable (optional) A data frame with selected bound variabl
#' @param min_bound (optional) A minimum bound value that applies to all clusters
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL "automatic regionalization with initial seed location"
#' @param initial_regions (optional) The initial regions that the local search starts with. Default is empty. means the local search starts with a random process to "grow" clusters
#' @param scale_method (optional) One of the scaling methods ('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust') to apply on input data. Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param random_seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage)
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares", "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' azp_clusters <- azp_tabu(5, queen_w, data, tabu_length=10, conv_tabu=10)
#' azp_clusters
#' }
#' @export
azp_tabu<- function(p, w, df, tabu_length=10, conv_tabu=10, bound_variable=data.frame(), min_bound=0, inits=0, initial_regions=vector('numeric'), scale_method="standardize", distance_method="euclidean", random_seed=123456789, rdist=numeric()) {
  if (p < 0) {
    stop("The p should be a positive integer number.")
  }
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }
  if (inherits(df, "data.frame") == FALSE) {
    stop("The input data needs to be a data.frame.")
  }

  n_vars <- length(df)

  if (inherits(df, "sf")) {
    n_vars <- n_vars - 1
  }

  if (n_vars < 1) {
    stop("The data.frame is empty.")
  }

  scale_methods <- c('raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust')
  if (!(scale_method %in% scale_methods)) {
    stop("The scale_method has to be one of {'raw', 'standardize', 'demean', 'mad', 'range_standardize', 'range_adjust'}")
  }

  if (distance_method != "euclidean" && distance_method != "manhattan") {
    stop("The distance method needs to be either 'euclidean' or 'manhattan'.")
  }

  bound_values <- vector('numeric')
  if (length(bound_variable) > 0)  {
    bound_values <- bound_variable[[1]]
  }

  return(p_azp_tabu(p, w$GetPointer(), df, n_vars, tabu_length, conv_tabu, bound_values, min_bound, inits, initial_regions, scale_method, distance_method, random_seed, rdist))
}

############################################################
#' @title Spatial Validation
#' @description Spatial validation provides a collection of validation measures including
#' 1. fragmentations (entropy, simpson), 2. join count ratio, 3. compactness (isoperimeter quotient)
#' and 4. diameter.
#' @param sf_obj An sf (simple feature) object
#' @param clusters A cluster classification variable (categorical values from a dataframe or values returned from cluster functions)
#' @param w An instance of Weight class
#' @return A list with names "Is Spatially Constrained", "Fragmentation", "Join Count Ratio",
#' "Compactness", and "Diameter".
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' clusters <- skater(5, queen_w, data)
#' results <- spatial_validation(guerry, clusters, queen_w)
#' results
#' }
#' @export
spatial_validation <- function(sf_obj, clusters, w) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }

  geoda_obj <- getGeoDaObj(sf_obj)

  if (geoda_obj$GetNumObs() <=0) {
    stop("gda object is not valid.")
  }

  return (p_spatialvalidation(geoda_obj$GetPointer(), clusters, w$GetPointer()));
}

############################################################
#' @title Join Count Ratio
#' @description Join count ratio is the join counts, the number of times a category is surrounded
#' by neighbors of the same category, over the total number of neighbors after converting
#' each category to a dummy variable.
#' @param clusters A cluster classification variable (categorical values from a dataframe or values returned from cluster functions)
#' @param w An instance of Weight class
#' @return A data.frame with names "Cluster", "N", "Neighbors", "Join Count", "Ratio"
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' queen_w <- queen_weights(guerry)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' clusters <- skater(5, queen_w, data)
#' results <- join_count_ratio(clusters, queen_w)
#' results
#' }
#' @export
join_count_ratio<- function(clusters, w) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }

  return (p_joincount_ratio(clusters, w$GetPointer()));
}

############################################################
#' @title Make Spatial
#' @description Make spatially constrained clusters from spatially non-constrained clusters
#' using the contiguity information from the input weights
#' @param clusters A cluster classification variable (categorical values from a dataframe or values returned from cluster functions)
#' @param w An instance of Weight class
#' @return A vector of categorical values (cluster classification)
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry <- st_read(guerry_path)
#' data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
#' clusters <- kmeans(5, data)
#' queen_w <- queen_weights(guerry)
#' results <- make_spatial(clusters, queen_w)
#' results
#' }
#' @export
make_spatial <- function(clusters, w) {
  if (w$num_obs < 1) {
    stop("The weights is not valid.")
  }

  if (w$num_obs != length(clusters)) {
    stop("The weights doesn not match with the size of input clusters.")
  }

  return (p_make_spatial(clusters, w$GetPointer()));
}
