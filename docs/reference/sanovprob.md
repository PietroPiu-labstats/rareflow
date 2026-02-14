# Sanov Probability Bound

Computes the classical Sanov upper bound:

## Usage

``` r
sanovprob(Q, P, n)
```

## Arguments

- Q:

  Observed empirical distribution.

- P:

  True distribution.

- n:

  Sample size.

## Value

A numeric upper bound.

## Details

\$\$P(Q_n \approx Q) \le \exp\\-n \\ KL(Q \\ P)\\\$\$

where Q is the empirical distribution and P is the true distribution.
