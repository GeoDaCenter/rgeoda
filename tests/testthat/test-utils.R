context("utils.R")


testthat::test_that('utils', {
    vec1 <- c(1,2,3,4,5,6,7)
    vec2 <- c(1.1,2.2,3.3,4.4,5.5)
    data <- list(vec1, vec2)

    testthat::expect_equal(
        demean(data),
        list(c(-3,-2, -1,  0,  1,  2, 3), c(-2.2, -1.1,  0.0,  1.1,  2.2))
    )
    testthat::expect_equal(
        natural_breaks(2, vec1),
        c(5)
    )
})

