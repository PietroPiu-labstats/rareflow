#' Radial Normalizing Flow (1D)
#'
#' A readable implementation of a 1-dimensional radial flow:
#'
#'   z_K = z_0 + beta / (alpha + |z_0 - z_ref|) * (z_0 - z_ref)
#'
#' where:
#'   - z_ref is a reference point
#'   - alpha > 0 ensures numerical stability
#'   - beta controls the strength of the radial deformation
#'
#' The log-determinant is computed analytically for the 1D case.
#'
#' @param z_ref Reference point for the radial transformation.
#' @param alpha Positive scalar controlling the denominator.
#' @param beta Scalar controlling the magnitude of the deformation.
#'
#' @return A flow model object with methods:
#'   - sampleq(n)
#'   - logq(z0)
#'   - applyflow(z0)
#' @export
radialflowmodel <- function(z_ref, alpha, beta) {

  # ------------------------------------------------------------
  # Forward transformation
  # ------------------------------------------------------------
  apply_flow <- function(z0) {
    r <- abs(z0 - z_ref)
    h <- beta / (alpha + r)
    zK <- z0 + h * (z0 - z_ref)
    zK
  }

  # ------------------------------------------------------------
  # Log-determinant of the Jacobian (1D)
  #
  #   d/dz [ z + h(z) * (z - z_ref) ]
  #   = 1 + h + (z - z_ref) * h'(z)
  #
  # where:
  #   h(z) = beta / (alpha + |z - z_ref|)
  #   h'(z) = -beta * sign(z - z_ref) / (alpha + |z - z_ref|)^2
  # ------------------------------------------------------------
  logdet <- function(z0) {
    r <- abs(z0 - z_ref)
    h <- beta / (alpha + r)
    hprime <- -beta / (alpha + r)^2
    deriv <- 1 + h + (z0 - z_ref) * sign(z0 - z_ref) * hprime
    log(abs(deriv))
  }

  # ------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------
  list(
    type = "radial",
    dim  = 1,
    z_ref = z_ref,
    alpha = alpha,
    beta  = beta,

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
