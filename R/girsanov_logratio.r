#' Girsanov Log-Ratio for Drift-Tilted Diffusions
#'
#' Computes the Radon–Nikodym derivative (log form) associated with a
#' Girsanov change of measure for an SDE:
#'
#'   \eqn{dX_t = b(X_t)\, dt + dW_t}
#'
#' tilted by an alternative drift:
#'
#'   \eqn{dX_t = (b(X_t) + \theta_t)\, dt + dW_t}.
#'
#' The log-likelihood ratio is:
#'
#'   \eqn{\log \frac{dQ}{dP}
#'     = \sum_t \left( \theta_t W_{t} - \frac{1}{2}\theta_t^2\, dt \right)}.
#'
#' This function returns the log-ratio for a given path of drift tilts
#' `theta_path`, Brownian increments `Winc`, and time step `dt`.
#'
#' @param theta_path Numeric vector of drift tilts θ_t.
#' @param Winc Numeric vector of Brownian increments ΔW_t.
#' @param dt Time step size.
#'
#' @return A numeric log-likelihood ratio.
#' @export
girsanov_logratio <- function(theta_path, Winc, dt) {

  # ------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------
  stopifnot(is.numeric(theta_path))
  stopifnot(is.numeric(Winc))
  stopifnot(length(theta_path) == length(Winc))
  stopifnot(dt > 0)

  # ------------------------------------------------------------
  # Girsanov log-likelihood ratio
  # ------------------------------------------------------------
  term1 <- sum(theta_path * Winc)
  term2 <- -0.5 * sum(theta_path^2) * dt

  term1 + term2
}
