#' Sigmoid Activation
#'
#' Standard logistic sigmoid function.
#'
#' @param x Numeric input.
#' @return Numeric value between 0 and 1.
#' @keywords internal
sigmoid <- function(x) 1 / (1 + exp(-x))


#' Neural Two-Separator Likelihood
#'
#' Constructs a simple 3-category likelihood model based on a latent vector `z`.
#' The likelihood is defined by two logistic separators:
#'
#'   p1 = sigmoid(mean(z) - a)
#'   p2 = sigmoid(mean(z) + a)
#'
#' producing a 3-class probability vector:
#'
#'   (1 - p1, p1 - p2, p2)
#'
#' This likelihood is useful for toy neural classification models or
#' simple latent-to-categorical mappings.
#'
#' @param a Separation parameter controlling the spacing between the two logits.
#'
#' @return A function mapping a latent vector `z` to a probability vector.
#' @export
makeneurolik <- function(a = 0.3) {
  function(z) {
    eta <- mean(z)
    p1 <- sigmoid(eta - a)
    p2 <- sigmoid(eta + a)
    out <- c(1 - p1, p1 - p2, p2)
    pmax(out / sum(out), 1e-12)
  }
}


#' Biological Two-Separator Likelihood
#'
#' Same structure as `makeneurolik()`, but with a different default separation
#' parameter. This can be used to model simple biological switching systems
#' or coarse-grained gene expression states.
#'
#' @param a Separation parameter controlling the spacing between the two logits.
#'
#' @return A function mapping a latent vector `z` to a probability vector.
#' @export
makebiolik <- function(a = 0.2) {
  function(z) {
    eta <- mean(z)
    p1 <- sigmoid(eta - a)
    p2 <- sigmoid(eta + a)
    out <- c(1 - p1, p1 - p2, p2)
    pmax(out / sum(out), 1e-12)
  }
}
