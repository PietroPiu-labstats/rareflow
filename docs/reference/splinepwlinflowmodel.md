# Piecewise-Linear Spline Flow (Monotone)

A readable implementation of a monotone piecewise-linear spline flow.
Each dimension is transformed independently using learned spline
parameters.

## Usage

``` r
splinepwlinflowmodel(d = 2, K = 8, theta = NULL)
```

## Arguments

- d:

  Dimension of the latent space.

- K:

  Number of spline bins.

- theta:

  Optional parameter vector. If NULL, random initialization.

## Value

A flow model object with methods:

- sampleq(n)

- logq(z0)

- applyflow(z0)

## Details

The spline flow uses:

- \\K\\ bins with learned widths \\w\\ and heights \\h\\

- a softmax transformation to ensure positivity and normalization

- a sigmoid reparameterization for numerical stability

The flow is invertible and differentiable almost everywhere.
