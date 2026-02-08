#' Freidlin–Wentzell Action Functional
#'
#' Computes the discrete Freidlin–Wentzell action for a path \eqn{\phi(t)}
#' represented as a matrix of size T × d. The continuous action is:
#'
#' \deqn{
#'   I[\phi] = \frac{1}{2} \int_0^T \| \dot{\phi}(t) - b(\phi(t)) \|^2 dt,
#' }
#'
#' and the discrete approximation is:
#'
#' \deqn{
#'   I \approx \frac{1}{2} \sum_{t=1}^{T-1}
#'     \| (\phi_{t+1} - \phi_t)/dt - b(\phi_t) \|^2 \, dt.
#' }
#'
#' @param phi Matrix of path values (T × d).
#' @param drift Drift function \eqn{b(x)} returning a numeric vector.
#' @param dt Time step.
#'
#' @return Numeric action value.
#' @export
Freidlin_Wentzell_action <- function(phi, drift, dt) {

  stopifnot(is.matrix(phi))
  stopifnot(is.function(drift))
  stopifnot(dt > 0)

  T <- nrow(phi)
  d <- ncol(phi)

  # (T-1) × d
  dphi <- apply(phi, 2, diff)
  dphi <- matrix(dphi, ncol = d)

  # (T-1) × d
  bphi <- matrix(NA_real_, nrow = T - 1, ncol = d)
  for (i in 1:(T - 1)) {
    out <- drift(phi[i, ])
    out <- as.numeric(out)
    if (length(out) != d) out <- rep(out, d)
    bphi[i, ] <- out
  }

  sum(rowSums((dphi / dt - bphi)^2)) * dt / 2
}


#' Freidlin–Wentzell Quasi-Potential via Path Minimization
#'
#' Computes an approximate Freidlin–Wentzell quasi-potential between two
#' points \eqn{x_0} and \eqn{x_1} by minimizing the FW action functional
#' over discretized paths.
#'
#' The algorithm:
#'   1. Initializes a straight-line path between x0 and x1.
#'   2. Performs simple gradient descent on the FW action.
#'
#' This is a naive but effective illustrative method for low-dimensional
#' systems. More advanced solvers (string method, MAM, etc.) can be plugged in.
#'
#' @param x0 Starting point (numeric vector).
#' @param x1 Target point (numeric vector).
#' @param drift Drift function b(x).
#' @param T Number of time steps.
#' @param dt Time step.
#' @param niter Number of gradient descent iterations.
#' @param stepsize Gradient descent step size.
#'
#' @return A list with:
#'   - path: matrix of size T × d
#'   - action: FW action of the optimized path
#'
#' @export
FW_quasipotential <- function(x0,
                              x1,
                              drift,
                              T = 200,
                              dt = 0.01,
                              niter = 200,
                              stepsize = 0.1) {

  stopifnot(is.numeric(x0), is.numeric(x1))
  stopifnot(length(x0) == length(x1))
  stopifnot(is.function(drift))
  stopifnot(T > 2, dt > 0)

  d <- length(x0)

  # Initialize path
  phi <- matrix(0, T, d)
  phi[1, ] <- x0
  phi[T, ] <- x1

  for (t in 2:(T - 1)) {
    phi[t, ] <- x0 + (t - 1) / (T - 1) * (x1 - x0)
  }

  # Gradient descent
  for (iter in 1:niter) {
    grad <- matrix(0, T, d)

    for (t in 2:(T - 1)) {
      dphi <- (phi[t + 1, ] - phi[t, ]) - (phi[t, ] - phi[t - 1, ])
      bval <- drift(phi[t, ])
      bval <- as.numeric(bval)
      if (length(bval) != d) bval <- rep(bval, d)
      grad[t, ] <- dphi - bval
    }

    phi <- phi - stepsize * grad
  }

  list(
    path = phi,
    action = Freidlin_Wentzell_action(phi, drift, dt)
  )
}
