#' Fit a Flow-Based Variational Posterior for Sanov Inference
#'
#' This function fits a variational posterior q(z) using a chosen normalizing flow.
#' The objective is the ELBO:
#'
#'   \eqn{\log p(x \mid z) + \log p(z) - \log q(z)}
#'
#' The flow parameters (theta) are optimized via `optim()` when applicable
#' (MAF and spline flows). Planar and radial flows have no trainable parameters.
#'
#' @param observed Observed empirical distribution Q (probability vector).
#' @param states Optional vector of category names.
#' @param flowtype One of "maf", "splinepwlin", "planar", "radial".
#' @param flowspec A list specifying structural parameters (d, K, etc.).
#' @param inittheta Optional initial parameter vector for trainable flows.
#' @param pxgivenz A function mapping latent z to a categorical pmf.
#' @param nmc Number of Monte Carlo samples for ELBO estimation.
#' @param control List of control parameters passed to `optim()`.
#'
#' @return A list containing:
#'   - flow: the fitted flow model
#'   - elbo: final ELBO value
#'   - theta: optimized parameter vector (if applicable)
#'   - convergence: optim() convergence code
#'
#' @details
#' This function performs generic variational inference using a chosen
#' normalizing flow and a user-provided likelihood `pxgivenz`.
#'
#' For specialized rare-event inference using:
#'   - Girsanov change of measure
#'   - Freidlin–Wentzell quasi-potential
#'
#' see the wrapper functions:
#'   - `fitflow_girsanov()`
#'   - `fitflow_FW()`
#'
#' These wrappers construct a tilted likelihood and then call
#' `fitflowvariational()` internally.
#'
#' @export
fitflowvariational <- function(observed,
                               states = NULL,
                               flowtype = c("maf", "splinepwlin", "planar", "radial"),
                               flowspec = list(),
                               inittheta = NULL,
                               pxgivenz,
                               nmc = 256,
                               control = list()) {

  # ------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------
  flowtype <- match.arg(flowtype)
  stopifnot(is.numeric(observed), abs(sum(observed) - 1) < 1e-8)
  stopifnot(is.function(pxgivenz))
  stopifnot(is.list(flowspec))

  # ------------------------------------------------------------
  # Initialize flow (with or without theta)
  # ------------------------------------------------------------
  flow0 <- switch(
    flowtype,

    planar = makeflow("planar", list(
      u = flowspec$u %||% 0.05,
      w = flowspec$w %||% 0.05,
      b = flowspec$b %||% 0
    )),

    radial = makeflow("radial", list(
      z_ref = flowspec$z_ref %||% 0,
      alpha = flowspec$alpha %||% 1,
      beta  = flowspec$beta  %||% 0.05
    )),

    maf = makeflow("maf", list(
      d = flowspec$d %||% 3,
      K = flowspec$K %||% 2,
      theta = inittheta
    )),

    splinepwlin = makeflow("splinepwlin", list(
      d = flowspec$d %||% 2,
      K = flowspec$K %||% 8,
      theta = inittheta
    ))
  )

  # ------------------------------------------------------------
  # Case 1: planar/radial → no trainable parameters
  # ------------------------------------------------------------
  if (flowtype %in% c("planar", "radial")) {
    elbo0 <- ELBOflow(flow0, observed, pxgivenz, nmc)
    return(list(
      flow = flow0,
      elbo = elbo0,
      theta = NULL,
      convergence = 0
    ))
  }

  # ------------------------------------------------------------
  # Case 2: trainable flows (MAF, spline)
  # ------------------------------------------------------------
  theta0 <- flow0$theta

  # Objective for optim(): negative ELBO
  objective <- function(theta) {
    flow <- makeflow(flowtype, list(
      d = flowspec$d,
      K = flowspec$K,
      theta = theta
    ))
    -ELBOflow(flow, observed, pxgivenz, nmc)
  }

  opt <- optim(
    par = theta0,
    fn = objective,
    method = "BFGS",
    control = modifyList(list(maxit = 300), control)
  )

  # Build final flow with optimized parameters
  flow_final <- makeflow(flowtype, list(
    d = flowspec$d,
    K = flowspec$K,
    theta = opt$par
  ))

  # ------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------
  list(
    flow = flow_final,
    elbo = -opt$value,
    theta = opt$par,
    convergence = opt$convergence
  )
}
