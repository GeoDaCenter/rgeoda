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

testthat::test_that("local_losh", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)

    # Test with a = 2 (default)
    losh2 <- local_losh(queen_w, guerry["Crm_prs"])
    hi_vals2 <- lisa_values(losh2)
    pvals2 <- lisa_pvalues(losh2)
    lbls <- lisa_labels(losh2)
    nn <- lisa_num_nbrs(losh2)

    # Values compared against spdep results
    testthat::expect_equal(hi_vals2[[1]], 0.4256315, tolerance = 1e-6)
    testthat::expect_equal(hi_vals2[[10]], 0.6160384, tolerance = 1e-6)
    testthat::expect_equal(lbls, c("Not significant", "Heterogeneous", "Homogeneous",
                                   "", "", "Undefined", "Isolated"))
    
    # Check neighbors for a few observations
    testthat::expect_equal(nn[[1]], 4)
    testthat::expect_equal(nn[[10]], 5)

    # Test with a = 1
    losh1 <- local_losh(queen_w, guerry["Crm_prs"], a = 1)
    hi_vals1 <- lisa_values(losh1)
    testthat::expect_equal(hi_vals1[[1]], 0.5942161, tolerance = 1e-6)

    # Test significance cutoff filter
    clusters_05 <- lisa_clusters(losh2, cutoff = 0.05)
    clusters_01 <- lisa_clusters(losh2, cutoff = 0.01)
    
    # Check that fewer observations are significant with tighter cutoff
    testthat::expect_true(sum(clusters_01 != 0) <= sum(clusters_05 != 0))
    
    # Check labels mapping
    # 1: Heterogeneous (Hi > 1), 2: Homogeneous (Hi < 1)
    for (i in 1:10) {
        if (clusters_05[[i]] == 1) {
            testthat::expect_true(hi_vals2[[i]] > 1.0)
        } else if (clusters_05[[i]] == 2) {
            testthat::expect_true(hi_vals2[[i]] < 1.0)
        }
    }

    # Test with kernel weights (which include self-links)
    kw <- kernel_knn_weights(guerry, k = 6, kernel_method = "gaussian")
    losh_kw <- local_losh(kw, guerry["Crm_prs"])
    kw_pvals <- lisa_pvalues(losh_kw)
    kw_hi <- lisa_values(losh_kw)
    testthat::expect_equal(length(kw_pvals), nrow(guerry))
    testthat::expect_true(all(!is.na(kw_pvals)))
    testthat::expect_true(all(!is.na(kw_hi)))
})

testthat::test_that("self-link invariance and isolate handling in local statistics", {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    queen_w <- queen_weights(guerry)
    n <- nrow(guerry)

    # Construct queen_w_self with self-links added to every observation
    queen_w_self <- create_weights(n)
    for (i in 1:n) {
        nbrs <- get_neighbors(queen_w, i)
        set_neighbors(queen_w_self, i, sort(unique(c(nbrs, i))))
    }
    update_weights(queen_w_self)

    # 1. Local Moran: complete permutation
    lm_w0 <- local_moran(queen_w, guerry["Crm_prs"], seed = 123456789, permutation_method = "complete")
    lm_w1 <- local_moran(queen_w_self, guerry["Crm_prs"], seed = 123456789, permutation_method = "complete")
    testthat::expect_equal(lisa_values(lm_w0), lisa_values(lm_w1))
    testthat::expect_equal(lisa_pvalues(lm_w0), lisa_pvalues(lm_w1))
    testthat::expect_equal(lisa_clusters(lm_w0), lisa_clusters(lm_w1))

    # 2. Local Moran: lookup permutation
    lm_w0_lookup <- local_moran(queen_w, guerry["Crm_prs"], seed = 123456789, permutation_method = "lookup")
    lm_w1_lookup <- local_moran(queen_w_self, guerry["Crm_prs"], seed = 123456789, permutation_method = "lookup")
    testthat::expect_equal(lisa_values(lm_w0_lookup), lisa_values(lm_w1_lookup))
    testthat::expect_equal(lisa_pvalues(lm_w0_lookup), lisa_pvalues(lm_w1_lookup))
    testthat::expect_equal(lisa_clusters(lm_w0_lookup), lisa_clusters(lm_w1_lookup))

    # 3. Local Join Count: complete permutation
    jc_w0 <- local_joincount(queen_w, guerry["TopCrm"], seed = 123456789, permutation_method = "complete")
    jc_w1 <- local_joincount(queen_w_self, guerry["TopCrm"], seed = 123456789, permutation_method = "complete")
    testthat::expect_equal(lisa_values(jc_w0), lisa_values(jc_w1))
    testthat::expect_equal(lisa_pvalues(jc_w0), lisa_pvalues(jc_w1))
    testthat::expect_equal(lisa_clusters(jc_w0), lisa_clusters(jc_w1))

    # 4. Local Join Count: lookup permutation
    jc_w0_lookup <- local_joincount(queen_w, guerry["TopCrm"], seed = 123456789, permutation_method = "lookup")
    jc_w1_lookup <- local_joincount(queen_w_self, guerry["TopCrm"], seed = 123456789, permutation_method = "lookup")
    testthat::expect_equal(lisa_values(jc_w0_lookup), lisa_values(jc_w1_lookup))
    testthat::expect_equal(lisa_pvalues(jc_w0_lookup), lisa_pvalues(jc_w1_lookup))
    testthat::expect_equal(lisa_clusters(jc_w0_lookup), lisa_clusters(jc_w1_lookup))

    # 5. Isolate handling: self-only neighbor (k=1 where neighbor is self)
    iso_w <- create_weights(6)
    set_neighbors(iso_w, 1, c(1))
    set_neighbors(iso_w, 2, c(3))
    set_neighbors(iso_w, 3, c(2, 4))
    set_neighbors(iso_w, 4, c(3, 5))
    set_neighbors(iso_w, 5, c(4, 6))
    set_neighbors(iso_w, 6, c(5))
    update_weights(iso_w)

    df_test <- data.frame(val = c(10.0, 1.0, 2.0, 3.0, 4.0, 5.0), bin = c(1, 0, 1, 0, 1, 0))
    lm_iso <- local_moran(iso_w, df_test["val"], seed = 123456789)
    testthat::expect_true(is.na(lisa_pvalues(lm_iso)[[1]]))
    testthat::expect_equal(lisa_clusters(lm_iso)[[1]], 6) # Isolated

    jc_iso <- local_joincount(iso_w, df_test["bin"], seed = 123456789)
    testthat::expect_true(is.na(lisa_pvalues(jc_iso)[[1]]))
    testthat::expect_equal(lisa_clusters(jc_iso)[[1]], 3) # Isolated

    losh_iso <- local_losh(iso_w, df_test["val"], seed = 123456789)
    testthat::expect_true(is.na(lisa_pvalues(losh_iso)[[1]]))
    testthat::expect_equal(lisa_clusters(losh_iso)[[1]], 6) # Isolated
})
