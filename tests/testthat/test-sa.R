context("sa.R")

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
