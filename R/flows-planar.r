#' Planar Normalizing Flow (1D)
#'
#' A simple and readable implementation of a 1-dimensional planar flow:
#'   z_K = z_0 + u * h(wz_0 + b)
#'
#' where:
#'   - h is a smooth activation (tanh)
#'   - the log-determinant is computed analytically
#'
#' This flow is mainly useful for pedagogical purposes or as a lightweight
#' building block in variational inference.
#'
#' @param u Scalar parameter controlling the magnitude of the deformation.
#' @param w Scalar parameter controlling the slope of the activation.
#' @param b Scalar bias term.
#'
#' @return A flow model object with methods:
#'   - sampleq(n)
#'   - logq(z0)
#'   - applyflow(z0)
#' @export
planarflowmodel <- function(u, w, b) {

  # ------------------------------------------------------------
  # Activation and derivative
  # ------------------------------------------------------------
  activate <- function(a) tanh(a)
  activate_prime <- function(a) 1 - tanh(a)^2

  # ------------------------------------------------------------
  # Forward transformation
  # ------------------------------------------------------------
  apply_flow <- function(z0) {
    a <- w * z0 + b
    zK <- z0 + u * activate(a)
    zK
  }

  # ------------------------------------------------------------
  # Log-determinant of the Jacobian
  # For 1D: log |1 + u * h'(w z + b) * w|
  # ------------------------------------------------------------
  logdet <- function(z0) {
    a <- w * z0 + b
    psi <- activate_prime(a) * w
    log(abs(1 + u * psi))
  }

  # ------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------
  list(
    type = "planar",
    dim  = 1,
    u = u,
    w = w,
    b = b,

    # sample from base Gaussian and transform
    sampleq = function(n = 1) {
      z0 <- rnorm(n)
      zK <- apply_flow(z0)
      list(
        z0 = matrix(z0, ncol = 1),
        zK = matrix(zK, ncol = 1)
      )
    },

    # compute log q(zK)
    logq = function(z0) {
      z <- as.numeric(z0[, 1])
      zK <- apply_flow(z)
      logq0 <- dnorm(z, log = TRUE)
      logq <- logq0 - logdet(z)
      list(
        zK = matrix(zK, ncol = 1),
        logq = logq
      )
    },

    # apply only the forward transformation
    applyflow = function(z0) {
      z <- as.numeric(z0[, 1])
      matrix(apply_flow(z), ncol = 1)
    }
  )
}
