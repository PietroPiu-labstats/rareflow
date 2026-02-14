# Fit Flow Using Freidlin–Wentzell Quasi-Potential

Computes the Freidlin–Wentzell quasi-potential between x0 and x1,
constructs a tilted likelihood proportional to exp(-V/eps), and fits a
flow-based variational posterior.

## Usage

``` r
fitflow_FW(
  observed,
  states = NULL,
  flowtype = "maf",
  flowspec = list(),
  inittheta = NULL,
  drift,
  x0,
  x1,
  T = 200,
  dt = 0.01,
  eps = 0.1,
  nmc = 256,
  control = list()
)
```

## Arguments

- observed:

  Empirical distribution Q.

- states:

  Optional category names.

- flowtype:

  Flow type.

- flowspec:

  Structural parameters for the flow.

- inittheta:

  Optional initial theta.

- drift:

  Drift function b(x).

- x0:

  Starting point.

- x1:

  Target point.

- T:

  Number of time steps.

- dt:

  Time step.

- eps:

  Noise strength (small parameter).

- nmc:

  Monte Carlo samples.

- control:

  Control list for optim().

## Value

Output of
[`fitflowvariational()`](https://pietropiu-labstats.github.io/rareflow/reference/fitflowvariational.md).

## Details

This is useful for rare-event inference in small-noise diffusions, where
the quasi-potential acts as an effective energy landscape.
