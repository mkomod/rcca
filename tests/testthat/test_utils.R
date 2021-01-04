testthat::test_that("Testing softmax", {
    testthat::expect_equal(
	rcca::soft_threshold(-3:3, 1.5),
	matrix(c(-1.5, -0.5, 0, 0, 0, 0.5, 1.5), ncol=1)
    )
})


testthat::test_that("Binary search", {
    # Binary search finds some threshold d such that the l1 norm of the
    # normalised softmax of vector w is equal to l
    w <- -3:3
    l <- runif(1, 1.5, 2.5)
    d <- rcca::binary_search(w, l, 1000, 1e-6)
    su <- rcca::soft_threshold(w, d)
    v <- su / norm(su, "F")
    l1.v <- sum(abs(v))
    testthat::expect_true(l1.v - l < 1e-6)
})
