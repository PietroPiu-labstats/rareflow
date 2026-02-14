# Evidence Lower Bound (ELBO) for Flow-Based Variational Inference

Computes the ELBO:

## Usage

``` r
ELBOflow(flow, Qobs, pxgivenz, nmc = 256)
```

## Arguments

- flow:

  A flow model created by
  [`makeflow()`](https://pietropiu-labstats.github.io/rareflow/reference/makeflow.md).

- Qobs:

  Observed empirical distribution (probability vector).

- pxgivenz:

  A function mapping a latent vector z to a categorical pmf.

- nmc:

  Number of Monte Carlo samples.

## Value

A numeric ELBO value.

## Details

\\\log p(x \\ z) + \log p(z) - \log q(z)\\

where:

- q(z) is the flow-based variational posterior

- p(z) is a standard Gaussian prior

- p(x \| z) is provided by the user as `pxgivenz`

The expectation is approximated via Monte Carlo sampling.

## Examples

``` r
f <- makeflow("planar", list(u = 0.1, w = 0.2, b = 0))
px <- function(z) c(0.3, 0.4, 0.3)
ELBOflow(f, Qobs = c(0.2, 0.5, 0.3), pxgivenz = px, nmc = 100)
#> [1] -1.061949
```
