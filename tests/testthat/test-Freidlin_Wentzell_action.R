test_that("Freidlin_Wentzell_action matches known value in 1D", {
  b <- function(x) c(as.numeric(x) - as.numeric(x)^3)
  phi <- matrix(seq(-1, 1, length.out = 200), ncol = 1)

  val <- Freidlin_Wentzell_action(phi, b, dt = 0.01)

  expect_equal(val, 1.080835, tolerance = 1e-6)
})
