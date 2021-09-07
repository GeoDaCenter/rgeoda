context("lisa.R")

testthat::test_that("local_moran", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    lm <- local_moran(queen_w, guerry["Crm_prs"])
    pvals <- lisa_pvalues(lm)
    lbls <- lisa_labels(lm)

    testthat::expect_equal(pvals[[1]], 0.197)
    testthat::expect_equal(lbls, c("Not significant", "High-High", "Low-Low",
                                   "Low-High", "High-Low", "Undefined",
                                   "Isolated"))
})

testthat::test_that("local_bimoran", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    lm <- local_bimoran(queen_w, guerry[c("Crm_prs", "Litercy")])
    lvals <- lisa_values(lm)
    pvals <- lisa_pvalues(lm)
    lbls <- lisa_labels(lm)

    testthat::expect_equal(lvals[[1]], 0.392663447638106)
    testthat::expect_equal(pvals[[1]], 0.2690)
    testthat::expect_equal(lbls, c("Not significant", "High-High", "Low-Low",
                                   "Low-High", "High-Low", "Undefined",
                                   "Isolated"))
})

testthat::test_that("local_multiquantilelisa", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    qsa <- local_multiquantilelisa(queen_w, guerry[c("Crm_prs", "Crm_prp")],
                                   k = c(4, 4), q = c(1, 1))
    pvals <- lisa_pvalues(qsa)

    testthat::expect_equal(pvals[[12]], 0.244)
})

testthat::test_that("localmoran_eb", {
    # NOTE: the data used for local moran eb statistics are meaningless,
    # just for testing
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    localeb <- local_moran_eb(queen_w, guerry[c("Crm_prs", "Pop1831")])

    pvals <- lisa_pvalues(localeb)

    testthat::expect_equal(pvals[[1]], 0.455)
})

testthat::test_that("neighbor_match_test", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    data <- guerry[c("Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants",
                     "Suicids")]
    nbr_test <- neighbor_match_test(data, 6)

    testthat::expect_equal(nbr_test["Probability"][[1]][[1]], 0.052638)
    testthat::expect_equal(nbr_test["Cardinality"][[1]][[1]], 2)
})

testthat::test_that("local_joincount", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    localjc_crm <- local_joincount(queen_w, guerry["TopCrm"])

    pvals <- lisa_pvalues(localjc_crm)

    testthat::expect_equal(pvals[[1]], 0.395)

})
