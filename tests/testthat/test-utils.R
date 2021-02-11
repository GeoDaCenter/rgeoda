context("utils.R")


testthat::test_that('natural_breaks', {
    guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
    guerry <- st_read(guerry_path)
    brks <- natural_breaks(5, guerry['Crm_prs'])
    testthat::expect_equal( brks, c(9474, 16722, 23316, 29872))
})

