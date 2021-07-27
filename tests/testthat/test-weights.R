context("weights.R")

testthat::test_that('read_gal', {
  guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
  guerry <- st_read(guerry_path)
  queen_w <- queen_weights(guerry)

  save_weights(queen_w, guerry['CODE_DE'], '/tmp/guerry1.gal')
  w <- read_gal('/tmp/guerry1.gal', guerry['CODE_DE'][[1]])

  testthat::expect_equal(queen_w$sparsity, w$sparsity)
})

testthat::test_that('read_gwt', {
  guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
  guerry <- st_read(guerry_path)
  knn_w <- knn_weights(guerry, k=4)

  save_weights(knn_w, guerry['CODE_DE'], '/tmp/guerry2.gwt')
  w <- read_gwt('/tmp/guerry2.gwt', guerry['CODE_DE'][[1]])

  testthat::expect_equal(knn_w$sparsity, w$sparsity)
})

testthat::test_that('read_swm', {
  w <- read_swm('/tmp/virginia_queen.swm')
  testthat::expect_equal(w$sparsity, 0.03168253)
})

testthat::test_that('as.matrix', {
  guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
  guerry <- st_read(guerry_path)

  library(spdep)

  w <- queen_weights(guerry)
  mw <- as.matrix(w)
  lw <- mat2listw(mw)

  knn_w <- knn_weights(guerry, k=4)
  knn_mw <- as.matrix(knn_w)
  knn_lw <- mat2listw(knn_mw)

})
