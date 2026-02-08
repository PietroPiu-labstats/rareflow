test_that("FW_quasipotential returns finite action and valid path", {
  b <- function(x) c(as.numeric(x) - as.numeric(x)^3)

  res <- FW_quasipotential(
    x0 = -1,
    x1 = 1,
    drift = b,
    T = 100,
    dt = 0.01,
    niter = 10,
    stepsize = 0.05
  )

  expect_true(is.list(res))
  expect_true(is.matrix(res$path))
  expect_true(is.finite(res$action))
  expect_equal(ncol(res$path), 1)
})
