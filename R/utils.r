#' Kullback–Leibler Divergence
#'
#' Computes the KL divergence D(Q || P) between two discrete probability
#' distributions Q and P.
#'
#' @param Q Observed distribution (probability vector).
#' @param P Baseline or reference distribution (probability vector).
#'
#' @return A numeric KL divergence value.
#' @export
KLdiv <- function(Q, P) {
  stopifnot(is.numeric(Q), is.numeric(P))
  stopifnot(length(Q) == length(P))
  idx <- Q > 0
  sum(Q[idx] * log(Q[idx] / P[idx]))
}


#' Sanov Probability Bound
#'
#' Computes the classical Sanov upper bound:
#'
#' @details
#' \deqn{P(Q_n \approx Q) \le \exp\{-n \, KL(Q \| P)\}}
#'
#' where Q is the empirical distribution and P is the true distribution.
#'
#' @param Q Observed empirical distribution.
#' @param P True distribution.
#' @param n Sample size.
#'
#' @return A numeric upper bound.
#' @export
sanovprob <- function(Q, P, n) {
  exp(-n * KLdiv(Q, P))
}


# -------------------------------------------------------------------
# Internal utilities
# -------------------------------------------------------------------

#' Sigmoid Activation
#'
#' Standard logistic sigmoid function.
#'
#' @param x Numeric input.
#' @return Numeric output in (0,1).
#' @keywords internal
sigmoid <- function(x) 1 / (1 + exp(-x))


#' Logit Transform
#'
#' Computes log(u / (1 - u)).
#'
#' @param u Numeric input in (0,1).
#' @return Numeric logit value.
#' @keywords internal
logit <- function(u) log(u) - log(1 - u)


#' Softmax Function
#'
#' Computes a numerically stable softmax.
#'
#' @param x Numeric vector.
#' @return Probability vector summing to 1.
#' @keywords internal
softmax <- function(x) {
  ex <- exp(x - max(x))
  ex / sum(ex)
}


#' Default Operator %||%
#'
#' Returns `b` if `a` is NULL, otherwise returns `a`.
#'
#' @param a Primary value.
#' @param b Default value.
#'
#' @return Either `a` or `b`.
#'
#' @name or_or
#' @aliases %||%
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
