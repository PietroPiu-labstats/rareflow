#' Piecewise-Linear Spline Flow (Monotone)
#'
#' A readable implementation of a monotone piecewise-linear spline flow.
#' Each dimension is transformed independently using:
#'
#'   - K bins with learned widths (w) and heights (h)
#'   - softmax ensures positivity and normalization
#'   - the transformation is applied in sigmoid space for stability
#'
#' The flow is invertible and differentiable almost everywhere.
#'
#' @param d Dimension of the latent space.
#' @param K Number of spline bins.
#' @param theta Optional parameter vector. If NULL, random initialization.
#'
#' @return A flow model object with methods:
#'   - sampleq(n)
#'   - logq(z0)
#'   - applyflow(z0)
#' @export
splinepwlinflowmodel <- function(d = 2, K = 8, theta = NULL) {

  # ------------------------------------------------------------
  # Parameter count:
  #   For each dimension j:
  #     - K widths (wl)
  #     - K heights (hl)
  #   Total = 2 * K * d
  # ------------------------------------------------------------
  if (is.null(theta)) {
    theta <- rnorm(2 * K * d, sd = 0.05)
  }

  # ------------------------------------------------------------
  # Unpack theta into per-dimension spline parameters
  # ------------------------------------------------------------
  unpack <- function(theta) {
    mats <- vector("list", d)
    off <- 0

    for (j in 1:d) {
      wl <- theta[(off + 1):(off + K)]
      off <- off + K

      hl <- theta[(off + 1):(off + K)]
      off <- off + K

      mats[[j]] <- list(wl = wl, hl = hl)
    }

    mats
  }

  mats <- unpack(theta)

  # ------------------------------------------------------------
  # Utilities: sigmoid and logit
  # ------------------------------------------------------------
  sigmoid <- function(x) 1 / (1 + exp(-x))
  logit   <- function(u) log(u) - log(1 - u)

  # ------------------------------------------------------------
  # Piecewise-linear spline mapping in [0,1]
  # ------------------------------------------------------------
  pwl_map <- function(u, wl, hl) {
    # normalize widths and heights
    w <- softmax(wl)
    h <- softmax(hl)

    # cumulative boundaries
    U <- c(0, cumsum(w))
    V <- c(0, cumsum(h))

    # clamp u to avoid numerical issues
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)

    # find bin index
    bin <- findInterval(u, U)
    bin[bin == length(U)] <- length(U) - 1
    i <- bin

    # bin boundaries
    uL <- U[i]
    uR <- U[i + 1]
    vL <- V[i]
    vR <- V[i + 1]

    slope <- (vR - vL) / (uR - uL)

    # forward mapping
    yu <- vL + slope * (u - uL)

    # log determinant (change of variables)
    logdet <- log(slope) +
      log(u) + log(1 - u) -
      log(yu) - log(1 - yu)

    list(yu = yu, logdet = logdet)
  }

  # ------------------------------------------------------------
  # Apply the spline flow to all dimensions
  # ------------------------------------------------------------
  apply_flow <- function(z0) {
    n <- nrow(z0)
    z <- matrix(0, n, d)
    logdet <- numeric(n)

    for (j in 1:d) {
      u <- sigmoid(z0[, j])
      step <- pwl_map(u, mats[[j]]$wl, mats[[j]]$hl)
      z[, j] <- logit(step$yu)
      logdet <- logdet + step$logdet
    }

    list(zK = z, logdet = logdet)
  }

  # ------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------
  list(
    type = "splinepwlin",
    dim  = d,
    K    = K,
    theta = theta,

    # sample from base Gaussian and transform
    sampleq = function(n = 1) {
      z0 <- matrix(rnorm(n * d), n, d)
      out <- apply_flow(z0)
      list(z0 = z0, zK = out$zK)
    },

    # compute log q(zK)
    logq = function(z0) {
      out <- apply_flow(z0)
      logq0 <- rowSums(dnorm(z0, log = TRUE))
      logq <- logq0 - out$logdet
      list(zK = out$zK, logq = logq)
    },

    # apply only the forward transformation
    applyflow = function(z0) apply_flow(z0)$zK
  )
}
