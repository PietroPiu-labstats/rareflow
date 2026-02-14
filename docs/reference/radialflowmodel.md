# Radial Normalizing Flow (1D)

A readable implementation of a 1-dimensional radial flow:

## Usage

``` r
radialflowmodel(z_ref, alpha, beta)
```

## Arguments

- z_ref:

  Reference point for the radial transformation.

- alpha:

  Positive scalar controlling the denominator.

- beta:

  Scalar controlling the magnitude of the deformation.

## Value

A flow model object with methods:

- sampleq(n)

- logq(z0)

- applyflow(z0)

## Details

z_K = z_0 + beta / (alpha + \|z_0 - z_ref\|) \* (z_0 - z_ref)

where:

- z_ref is a reference point

- alpha \> 0 ensures numerical stability

- beta controls the strength of the radial deformation

The log-determinant is computed analytically for the 1D case.
