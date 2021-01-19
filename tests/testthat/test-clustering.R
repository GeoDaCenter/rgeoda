context("clustering.R")

testthat::test_that('schc', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    clusters <- schc(5, queen_w, data, 'average')

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(clusters, data)
    ratio <- betweenss / totalss


    testthat::expect_equal( ratio, 0.21477112549551957699)

})

# NOTE!!!!!!!!!
# The results are computed using Boost library 1.58.0. To pass the following test cases
# , please install BH package version==1.58.0

testthat::test_that('azp_greedy', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    azp_clusters <- azp_greedy(5, queen_w, data)

    bound_vals <- guerry_df['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831

    azp_clusters <- azp_greedy(5, queen_w, data, bound_vals = bound_vals, min_bound = min_bound)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(azp_clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.37015555244290943)

})

testthat::test_that('azp_sa', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    azp_clusters <- azp_sa(5, queen_w, data, cooling_rate = 0.85, sa_maxit = 1)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(azp_clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.4302644093198012)

})

testthat::test_that('azp_tabu', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    azp_clusters <- azp_tabu(5, queen_w, data, tabu_length = 10, conv_tabu=10)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(azp_clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.42221739641363148)

})

testthat::test_that('maxp_greedy', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    bound_vals <- guerry_df['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831

    clusters <- maxp_greedy(queen_w, data, bound_vals, min_bound)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.42329309419590377)

})

testthat::test_that('maxp_sa', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    bound_vals <- guerry_df['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831

    clusters <- maxp_sa(queen_w, data, bound_vals, min_bound, cooling_rate = 0.85, sa_maxit = 1)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.46832569616103947)

})

testthat::test_that('maxp_tabu', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)
    data <- as.list(guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')])

    bound_vals <- guerry_df['Pop1831'][,1]
    min_bound <- 3236.67 # 10% of Pop1831

    clusters <- maxp_tabu(queen_w, data, bound_vals, min_bound, tabu_length = 10, conv_tabu=10)

    totalss <- total_sumofsquare( data )
    betweenss <- between_sumofsquare(clusters, data)
    ratio <- betweenss / totalss

    testthat::expect_equal( ratio, 0.4893668149272537)

})
