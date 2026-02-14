# Piecewise-Linear Spline Flow (Monotone)

A readable implementation of a monotone piecewise-linear spline flow.
Each dimension is transformed independently using:

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

- K bins with learned widths (w) and heights (h)

- softmax ensures positivity and normalization

- the transformation is applied in sigmoid space for stability

The flow is invertible and differentiable almost everywhere.
