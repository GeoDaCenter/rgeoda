context("clustering.R")

testthat::test_that("spatial_validataion", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]
    clusters <- skater(6, queen_w, data)

    validation <- spatial_validation(guerry, clusters[[1]], queen_w)

    testthat::expect_equal(validation$IsSpatiallyConstrained, TRUE)
    testthat::expect_equal(validation$Fragmentation$Entropy, 1.5302035777210896)
    testthat::expect_equal(validation$Fragmentation$`Entropy*`, 0.85402287751287753)
    testthat::expect_equal(validation$Fragmentation$Simpson, 0.25619377162629758)
    testthat::expect_equal(validation$Fragmentation$`Simpson*`, 1.5371626297577856)
    testthat::expect_equal(validation$JoinCountRatio$Ratio, c(0.8571429,0.8923077,0.5846154,0.5454545,0.3846154,0.6666667), tolerance=1e-6)
    testthat::expect_equal(validation$Compactness$IPC, c(0.009772352, 0.009914427, 0.029675045, 0.034800225, 0.046733291, 0.035828472))
    testthat::expect_equal(validation$Diameter$Ratio, c(0.2413793, 0.2500000, 0.36363636363, 0.3750000, 0.6000000, 0.5000000))
})

testthat::test_that("make_spatial", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    km_clusters <- c(5,2,5,1,3,6,2,5,3,1,5,3,4,5,4,4,4,4,2,6,5,1,3,1,3,3,4,1,1,1,2,1,6,4,1,1,2,4,1,5,5,1,3,1,1,1,2,2,3,2,2,2,2,4,3,4,2,2,2,3,5,1,5,1,3,3,3,2,2,2,3,3,3,3,4,2,1,1,1,1,6,6,4,2,3)
    clusters <- make_spatial(km_clusters, queen_w)

    testthat::expect_equal(clusters, c(1, 2, 5, 1, 1, 1,2,1,3,1,5,1,3,5,4,4,4,4,2,6,5,1,2,1,3,3,4,1,1,1,1,1,6,4,4,1,2,1,4,5,5,4,3,1,1,1,4,4,3,2,4,2,2,4,2,4,2,2,4,2,5,1,1,1,2,2,1,2,2,4,3,3,3,3,4,2,1,1,1,1,4,4,4,2,3))
})

testthat::test_that("skater", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]
    data$geometry = NULL
    df <- scale(data)
    dv <- as.vector(dist(df))
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]
    clusters <- skater(5, queen_w, data, rdist = dv)

    testthat::expect_equal(clusters[[5]], 0.3763086809)
})

testthat::test_that("schc", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    clusters <- schc(5, queen_w, data, "average")

    testthat::expect_equal(clusters[[5]], 0.2147711255)
})

# NOTE
# The previous results are computed using Boost library 1.58.0.
# The new results are computed using Boost library 1.81.0.1
# The differences are caused by the different implementation of
# boost::unordered_map: he keys in boost::unordered_map are not ordered and
# have different orders in the two Boost versions. This involves a different
# mechanism of randomness in max-p algorithm when picking which area or region
# to process.

testthat::test_that("azp_greedy", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    azp_clusters <- azp_greedy(5, queen_w, data)

    testthat::expect_equal(azp_clusters[[5]], 0.36, tolerance = 1e-3)

    bound_variable <- guerry["Pop1831"]
    min_bound <- 3236.67 # 10% of Pop1831

    azp_clusters <- azp_greedy(5, queen_w, data,
                               bound_variable = bound_variable,
                               min_bound = min_bound)

    testthat::expect_equal(azp_clusters[[5]], 0.417, tolerance = 1e-3)

})

testthat::test_that("azp_sa", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    azp_clusters <- azp_sa(5, queen_w, data, cooling_rate = 0.85, sa_maxit = 1)

    testthat::expect_equal(azp_clusters[[5]], 0.359, tolerance = 1e-3)
})

testthat::test_that("azp_tabu", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    azp_clusters <- azp_tabu(5, queen_w, data, tabu_length = 10, conv_tabu=10)

    testthat::expect_equal(azp_clusters[[5]], 0.4222174)

})

testthat::test_that("maxp_greedy", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    bound_vals <- guerry["Pop1831"]
    min_bound <- 3236.67 # 10% of Pop1831

    clusters <- maxp_greedy(queen_w, data, bound_vals, min_bound)

    testthat::expect_equal(clusters[[5]], 0.484, tolerance = 1e-3)
})

testthat::test_that("maxp_sa", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    bound_vals <- guerry["Pop1831"]
    min_bound <- 3236.67 # 10% of Pop1831

    clusters <- maxp_sa(queen_w, data, bound_vals, min_bound,
                        cooling_rate = 0.85, sa_maxit = 1)

    testthat::expect_equal(clusters[[5]], 0.496, tolerance = 1e-3)
})

testthat::test_that("maxp_tabu", {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]

    bound_vals <- guerry["Pop1831"]
    min_bound <- 3236.67 # 10% of Pop1831


    clusters <- maxp_tabu(queen_w, data, bound_vals, min_bound,
                          tabu_length = 10, conv_tabu = 10)

    testthat::expect_equal(clusters[[5]], 0.478, tolerance = 1e-3)

})
