context("clustering.R")

testthat::test_that('skater', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
    clusters <- skater(5, queen_w, data)

    testthat::expect_equal(clusters[[5]], 0.3763086809)
})

testthat::test_that('schc', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    clusters <- schc(5, queen_w, data, 'average')

    testthat::expect_equal(clusters[[5]], 0.2147711255)
})

# NOTE!!!!!!!!!
# The results are computed using Boost library 1.58.0. To pass the following test cases
# , please install BH package version==1.58.0

testthat::test_that('azp_greedy', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    azp_clusters <- azp_greedy(5, queen_w, data)

    testthat::expect_equal(azp_clusters[[5]], 0.3598541)

    bound_variable <- guerry['Pop1831']
    min_bound <- 3236.67 # 10% of Pop1831

    azp_clusters <- azp_greedy(5, queen_w, data, bound_variable = bound_variable, min_bound = min_bound)

    testthat::expect_equal(azp_clusters[[5]], 0.3980921835)

})

testthat::test_that('azp_sa', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    azp_clusters <- azp_sa(5, queen_w, data, cooling_rate = 0.85, sa_maxit = 1)

    testthat::expect_equal(azp_clusters[[5]], 0.4211363)
})

testthat::test_that('azp_tabu', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    azp_clusters <- azp_tabu(5, queen_w, data, tabu_length = 10, conv_tabu=10)

    testthat::expect_equal(azp_clusters[[5]], 0.4222174)

})

testthat::test_that('maxp_greedy', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    bound_vals <- guerry['Pop1831']
    min_bound <- 3236.67 # 10% of Pop1831

    #clusters <- maxp_greedy(queen_w, data, bound_vals, min_bound)

    #testthat::expect_equal(clusters[[5]], 0.4499671068)
})

testthat::test_that('maxp_sa', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    bound_vals <- guerry['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831

    #clusters <- maxp_sa(queen_w, data, bound_vals, min_bound, cooling_rate = 0.85, sa_maxit = 1)

    #testthat::expect_equal(clusters[[5]], 0.4585352223)
})

testthat::test_that('maxp_tabu', {
    library(sf)
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]

    bound_vals <- guerry['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831


    #clusters <- maxp_tabu(queen_w, data, bound_vals, min_bound, tabu_length = 10, conv_tabu=10)

    #testthat::expect_equal(clusters[[5]], 0.4893668149)

})
