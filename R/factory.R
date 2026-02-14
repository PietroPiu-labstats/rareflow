#' Flow Factory
#'
#' A unified constructor for all flow models in the rareflow package.
#' This function dispatches to the appropriate flow implementation
#' based on the `flowtype` argument.
#'
#' Supported flow types:
#'   - "planar"       -> planarflowmodel()
#'   - "radial"       -> radialflowmodel()
#'   - "maf"          -> mafflowmodel()
#'   - "splinepwlin"  -> splinepwlinflowmodel()
#'
#' @param flowtype Character string specifying the flow type.
#' @param params A named list of parameters required by the chosen flow.
#'
#' @return A flow model object with methods:
#'   - sampleq(n)
#'   - logq(z0)
#'   - applyflow(z0)
#'
#' @examples
#' # Create a planar flow
#' f <- makeflow("planar", list(u = 0.1, w = 0.2, b = 0))
#' s <- f$sampleq(10)
#'
#' # Create a 2D spline flow
#' f2 <- makeflow("splinepwlin", list(d = 2, K = 8))
#'
#' @export
makeflow <- function(flowtype,
                     params = list()) {

  flowtype <- match.arg(flowtype,
                        choices = c("planar", "radial", "maf", "splinepwlin"))

  # ------------------------------------------------------------
  # Dispatch to the appropriate flow constructor
  # ------------------------------------------------------------
  flow <- switch(
    flowtype,

    planar = {
      stopifnot(all(c("u", "w", "b") %in% names(params)))
      planarflowmodel(
        u = params$u,
        w = params$w,
        b = params$b
      )
    },

    radial = {
      stopifnot(all(c("z_ref", "alpha", "beta") %in% names(params)))
      radialflowmodel(
        z_ref = params$z_ref,
        alpha = params$alpha,
        beta  = params$beta
      )
    },

    maf = {
      # defaults for MAF
      d <- params$d %||% 3
      K <- params$K %||% 2
      theta <- params$theta %||% NULL

      mafflowmodel(
        d = d,
        K = K,
        theta = theta
      )
    },

    splinepwlin = {
      # defaults for spline flow
      d <- params$d %||% 2
      K <- params$K %||% 8
      theta <- params$theta %||% NULL

      splinepwlinflowmodel(
        d = d,
        K = K,
        theta = theta
      )
    }
  )

  flow
}
