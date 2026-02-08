#' Masked Autoregressive Flow (MAF)
#'
#' A readable and structured implementation of a Masked Autoregressive Flow.
#' This version supports:
#'   - arbitrary dimension d
#'   - K sequential flow steps
#'   - a single parameter vector theta containing all weights
#'
#' @param d Dimension of the latent space.
#' @param K Number of flow steps.
#' @param theta Optional parameter vector. If NULL, random initialization.
#'
#' @return A flow model object with methods:
#'   - sampleq(n)
#'   - logq(z0)
#'   - applyflow(z0)
#' @export
mafflowmodel <- function(d = 3, K = 2, theta = NULL) {

  # ------------------------------------------------------------
  # Helper: number of strictly lower-triangular entries in d x d
  # ------------------------------------------------------------
  lowercount <- function(d) d * (d - 1) / 2

  # ------------------------------------------------------------
  # Total number of parameters per flow step:
  #   - bs: d biases for scale network
  #   - Ls: lower triangular weights for scale network
  #   - bt: d biases for translation network
  #   - Lt: lower triangular weights for translation network
  # ------------------------------------------------------------
  perstep <- 2 * d + 2 * lowercount(d)

  # ------------------------------------------------------------
  # Initialize theta if not provided
  # ------------------------------------------------------------
  if (is.null(theta)) {
    theta <- rnorm(K * perstep, sd = 0.05)
  }

  # ------------------------------------------------------------
  # Unpack theta into K flow steps
  # Each step contains:
  #   bs (d), Ls (lower triangle), bt (d), Lt (lower triangle)
  # ------------------------------------------------------------
  unpack <- function(theta) {
    steps <- vector("list", K)
    off <- 0

    for (k in 1:K) {
      # scale biases
      bs <- theta[(off + 1):(off + d)]
      off <- off + d

      # scale lower-triangular weights
      Ls_vec <- theta[(off + 1):(off + lowercount(d))]
      off <- off + lowercount(d)

      # translation biases
      bt <- theta[(off + 1):(off + d)]
      off <- off + d

      # translation lower-triangular weights
      Lt_vec <- theta[(off + 1):(off + lowercount(d))]
      off <- off + lowercount(d)

      # build full matrices
      Ls <- matrix(0, d, d)
      Lt <- matrix(0, d, d)

      idx <- 1
      for (i in 2:d) {
        for (j in 1:(i - 1)) {
          Ls[i, j] <- Ls_vec[idx]
          Lt[i, j] <- Lt_vec[idx]
          idx <- idx + 1
        }
      }

      steps[[k]] <- list(bs = bs, Ls = Ls, bt = bt, Lt = Lt)
    }

    steps
  }

  steps <- unpack(theta)

  # ------------------------------------------------------------
  # Forward transformation for a single flow step
  # ------------------------------------------------------------
  stepforward <- function(x, step) {
    n <- nrow(x)
    y <- x
    logdet <- numeric(n)

    for (i in 1:d) {
      # previous dimensions (autoregressive)
      if (i == 1) {
        m <- matrix(0, n, 0)
      } else {
        m <- y[, 1:(i - 1), drop = FALSE]
      }

      # scale network: s_i = tanh(b_s + m %*% W_s)
      si <- step$bs[i]
      if (i > 1) si <- si + as.numeric(m %*% step$Ls[i, 1:(i - 1)])
      si <- tanh(si)

      # translation network: t_i = b_t + m %*% W_t
      ti <- step$bt[i]
      if (i > 1) ti <- ti + as.numeric(m %*% step$Lt[i, 1:(i - 1)])

      # affine transformation
      y[, i] <- exp(si) * y[, i] + ti

      # log determinant contribution
      logdet <- logdet + si
    }

    list(y = y, logdet = logdet)
  }

  # ------------------------------------------------------------
  # Apply all K flow steps
  # ------------------------------------------------------------
  applyflow <- function(z0) {
    z <- z0
    total_logdet <- rep(0, nrow(z))

    for (k in 1:K) {
      out <- stepforward(z, steps[[k]])
      z <- out$y
      total_logdet <- total_logdet + out$logdet
    }

    list(zK = z, logdet = total_logdet)
  }

  # ------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------
  list(
    type = "maf",
    dim = d,
    K = K,
    theta = theta,

    # sample from base Gaussian and transform
    sampleq = function(n = 1) {
      z0 <- matrix(rnorm(n * d), n, d)
      out <- applyflow(z0)
      list(z0 = z0, zK = out$zK)
    },

    # compute log q(zK)
    logq = function(z0) {
      out <- applyflow(z0)
      logq0 <- rowSums(dnorm(z0, log = TRUE))
      logq <- logq0 - out$logdet
      list(zK = out$zK, logq = logq)
    },

    # apply only the forward transformation
    applyflow = function(z0) applyflow(z0)$zK
  )
}
