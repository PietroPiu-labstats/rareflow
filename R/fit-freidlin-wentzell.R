#' Fit Flow Using Freidlin–Wentzell Quasi-Potential
#'
#' Computes the Freidlin–Wentzell quasi-potential between x0 and x1,
#' constructs a tilted likelihood proportional to exp(-V/eps),
#' and fits a flow-based variational posterior.
#'
#' This is useful for rare-event inference in small-noise diffusions,
#' where the quasi-potential acts as an effective energy landscape.
#'
#' @param observed Empirical distribution Q.
#' @param states Optional category names.
#' @param flowtype Flow type.
#' @param flowspec Structural parameters for the flow.
#' @param inittheta Optional initial theta.
#' @param drift Drift function b(x).
#' @param x0 Starting point.
#' @param x1 Target point.
#' @param T Number of time steps.
#' @param dt Time step.
#' @param eps Noise strength (small parameter).
#' @param nmc Monte Carlo samples.
#' @param control Control list for optim().
#'
#' @return Output of `fitflowvariational()`.
#' @export
fitflow_FW <- function(observed,
                       states = NULL,
                       flowtype = "maf",
                       flowspec = list(),
                       inittheta = NULL,
                       drift,
                       x0,
                       x1,
                       T = 200,
                       dt = 0.01,
                       eps = 0.1,
                       nmc = 256,
                       control = list()) {

  qp <- FW_quasipotential(x0, x1, drift, T = T, dt = dt)
  V <- qp$action

  pxgivenz <- function(z) {
    w <- exp(-V / eps)
    p <- rep(1, length(observed)) * w
    p / sum(p)
  }

  fitflowvariational(
    observed = observed,
    states = states,
    flowtype = flowtype,
    flowspec = flowspec,
    inittheta = inittheta,
    pxgivenz = pxgivenz,
    nmc = nmc,
    control = control
  )
}
