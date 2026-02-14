# Fit Flow with Girsanov-Tilted Likelihood

Applies a Girsanov change of measure to tilt the likelihood and then
fits a flow-based variational posterior using
[`fitflowvariational()`](https://pietropiu-labstats.github.io/rareflow/reference/fitflowvariational.md).

## Usage

``` r
fitflow_girsanov(
  observed,
  states = NULL,
  flowtype = "maf",
  flowspec = list(),
  inittheta = NULL,
  base_pxgivenz,
  theta_path,
  Winc,
  dt,
  nmc = 256,
  control = list()
)
```

## Arguments

- observed:

  Empirical distribution Q (probability vector).

- states:

  Optional category names.

- flowtype:

  Flow type ("maf", "splinepwlin", "planar", "radial").

- flowspec:

  Structural parameters for the flow.

- inittheta:

  Optional initial theta for trainable flows.

- base_pxgivenz:

  Likelihood \\p(x \mid z)\\ before tilting.

- theta_path:

  Drift-tilting function or vector for Girsanov.

- Winc:

  Brownian increments.

- dt:

  Time step.

- nmc:

  Monte Carlo samples.

- control:

  Control list for [`optim()`](https://rdrr.io/r/stats/optim.html).

## Value

Output of
[`fitflowvariational()`](https://pietropiu-labstats.github.io/rareflow/reference/fitflowvariational.md).

## Details

This is useful when the target distribution arises from a drift-tilted
diffusion process, where the Radon-Nikodym derivative is given by the
Girsanov theorem.
