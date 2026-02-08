#' Evidence Lower Bound (ELBO) for Flow-Based Variational Inference
#'
#' Computes the ELBO:
#'
#'   \eqn{\log p(x \mid z) + \log p(z) - \log q(z)}
#'
#' where:
#'   - q(z) is the flow-based variational posterior
#'   - p(z) is a standard Gaussian prior
#'   - p(x | z) is provided by the user as `pxgivenz`
#'
#' The expectation is approximated via Monte Carlo sampling.
#'
#' @param flow A flow model created by `makeflow()`.
#' @param Qobs Observed empirical distribution (probability vector).
#' @param pxgivenz A function mapping a latent vector z to a categorical pmf.
#' @param nmc Number of Monte Carlo samples.
#'
#' @return A numeric ELBO value.
#'
#' @examples
#' f <- makeflow("planar", list(u = 0.1, w = 0.2, b = 0))
#' px <- function(z) c(0.3, 0.4, 0.3)
#' ELBOflow(f, Qobs = c(0.2, 0.5, 0.3), pxgivenz = px, nmc = 100)
#'
#' @export
ELBOflow <- function(flow,
                     Qobs,
                     pxgivenz,
                     nmc = 256) {

  # ------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------
  stopifnot(is.list(flow))
  stopifnot(is.numeric(Qobs), abs(sum(Qobs) - 1) < 1e-8)
  stopifnot(is.function(pxgivenz))
  stopifnot(nmc > 0)

  # ------------------------------------------------------------
  # Sample from q(z0) and transform through the flow
  # ------------------------------------------------------------
  s <- flow$sampleq(nmc)
  z0 <- s$z0

  q <- flow$logq(z0)
  zK <- q$zK
  logq <- q$logq

  # ------------------------------------------------------------
  # Compute log-likelihood term: E_q[ log p(x | z) ]
  # ------------------------------------------------------------
  loglik <- apply(zK, 1, function(z) {
    px <- pxgivenz(as.numeric(z))
    px <- pmax(px, 1e-12)  # avoid log(0)
    sum(Qobs * log(px))
  })

  # ------------------------------------------------------------
  # Prior term: log p(z) for standard Gaussian
  # ------------------------------------------------------------
  logpz <- rowSums(dnorm(zK, log = TRUE))

  # ------------------------------------------------------------
  # Monte Carlo estimate of ELBO
  # ------------------------------------------------------------
  mean(loglik + logpz - logq)
}
