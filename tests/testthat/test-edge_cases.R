test_that("Freidlin_Wentzell_action handles minimal T", {
  b <- function(x) x
  phi <- matrix(c(0, 1), ncol = 1)
  expect_silent(Freidlin_Wentzell_action(phi, b, dt = 0.1))
})
