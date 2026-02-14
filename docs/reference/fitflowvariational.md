# Fit a Flow-Based Variational Posterior for Sanov Inference

This function fits a variational posterior q(z) using a chosen
normalizing flow. The objective is the ELBO:

## Usage

``` r
fitflowvariational(
  observed,
  states = NULL,
  flowtype = c("maf", "splinepwlin", "planar", "radial"),
  flowspec = list(),
  inittheta = NULL,
  pxgivenz,
  nmc = 256,
  control = list()
)
```

## Arguments

- observed:

  Observed empirical distribution Q (probability vector).

- states:

  Optional vector of category names.

- flowtype:

  One of "maf", "splinepwlin", "planar", "radial".

- flowspec:

  A list specifying structural parameters (d, K, etc.).

- inittheta:

  Optional initial parameter vector for trainable flows.

- pxgivenz:

  A function mapping latent z to a categorical pmf.

- nmc:

  Number of Monte Carlo samples for ELBO estimation.

- control:

  List of control parameters passed to
  [`optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list containing:

- flow: the fitted flow model

- elbo: final ELBO value

- theta: optimized parameter vector (if applicable)

- convergence: optim() convergence code

## Details

\\\log p(x \| z) + \log p(z) - \log q(z)\\

The flow parameters (theta) are optimized via
[`optim()`](https://rdrr.io/r/stats/optim.html) when applicable (MAF and
spline flows). Planar and radial flows have no trainable parameters.

This function performs generic variational inference using a chosen
normalizing flow and a user-provided likelihood `pxgivenz`.

For specialized rare-event inference using:

- Girsanov change of measure

- Freidlin–Wentzell quasi-potential

see the wrapper functions:

- [`fitflow_girsanov()`](https://pietropiu-labstats.github.io/rareflow/reference/fitflow_girsanov.md)

- [`fitflow_FW()`](https://pietropiu-labstats.github.io/rareflow/reference/fitflow_FW.md)

These wrappers construct a tilted likelihood and then call
`fitflowvariational()` internally.
