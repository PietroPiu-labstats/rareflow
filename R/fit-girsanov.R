#' Fit Flow with Girsanov-Tilted Likelihood
#'
#' Applies a Girsanov change of measure to tilt the likelihood and then
#' fits a flow-based variational posterior using `fitflowvariational()`.
#'
#' This is useful when the target distribution arises from a drift-tilted
#' diffusion process, where the Radon–Nikodym derivative is given by the
#' Girsanov theorem.
#'
#' @param observed Empirical distribution Q (probability vector).
#' @param states Optional category names.
#' @param flowtype Flow type ("maf", "splinepwlin", "planar", "radial").
#' @param flowspec Structural parameters for the flow.
#' @param inittheta Optional initial theta for trainable flows.
#' @param base_pxgivenz Likelihood p(x|z) before tilting.
#' @param theta_path Drift-tilting function or vector for Girsanov.
#' @param Winc Brownian increments.
#' @param dt Time step.
#' @param nmc Monte Carlo samples.
#' @param control Control list for optim().
#'
#' @return Output of `fitflowvariational()`.
#' @export
fitflow_girsanov <- function(observed,
                             states = NULL,
                             flowtype = "maf",
                             flowspec = list(),
                             inittheta = NULL,
                             base_pxgivenz,
                             theta_path,
                             Winc,
                             dt,
                             nmc = 256,
                             control = list()) {

  # Compute Girsanov log-weight
  logw <- girsanov_logratio(theta_path, Winc, dt)
  w <- exp(logw)

  # Tilted likelihood
  px_tilted <- function(z) {
    px <- base_pxgivenz(z)
    px <- pmax(px, 1e-12)
    p <- px * w
    p / sum(p)
  }

  fitflowvariational(
    observed = observed,
    states = states,
    flowtype = flowtype,
    flowspec = flowspec,
    inittheta = inittheta,
    pxgivenz = px_tilted,
    nmc = nmc,
    control = control
  )
}
