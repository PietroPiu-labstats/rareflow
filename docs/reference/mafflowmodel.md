# Masked Autoregressive Flow (MAF)

A readable and structured implementation of a Masked Autoregressive
Flow. This version supports:

- arbitrary dimension d

- K sequential flow steps

- a single parameter vector theta containing all weights

## Usage

``` r
mafflowmodel(d = 3, K = 2, theta = NULL)
```

## Arguments

- d:

  Dimension of the latent space.

- K:

  Number of flow steps.

- theta:

  Optional parameter vector. If NULL, random initialization.

## Value

A flow model object with methods:

- sampleq(n)

- logq(z0)

- applyflow(z0)
