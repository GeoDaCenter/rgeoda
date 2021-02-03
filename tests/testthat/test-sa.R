context("lisa.R")

testthat::test_that('localmoran_eb', {
    # NOTE: the data used for local moran eb statistics are meaningless, just for testing
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    queen_w <- queen_weights(guerry)
    guerry_df <- as.data.frame(guerry)
    crm <- guerry_df['Crm_prs'][,1]
    pop <- guerry_df['Pop1831'][,1]

    localeb <- local_moran_eb(queen_w, crm, pop)

    pvals <- lisa_pvalues(localeb)

    testthat::expect_equal( pvals[[1]], 0.455)
})

testthat::test_that('neighbor_match_test', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    data <- guerry_df[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
    nbr_test <- neighbor_match_test(guerry, 6, data)

    testthat::expect_equal( nbr_test['Probability'][[1]][[1]], 0.052638)
    testthat::expect_equal( nbr_test['Cardinality'][[1]][[1]], 2)
})

testthat::test_that('local_joincount', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- geoda_open(guerry_path)
    guerry_df <- as.data.frame(guerry)
    queen_w <- queen_weights(guerry)

    top_crm <- guerry_df['TopCrm'][,1]
    localjc_crm <- local_joincount(queen_w, top_crm)

    pvals <- lisa_pvalues(localjc_crm)

    testthat::expect_equal( pvals[[1]], 0.395)

})
