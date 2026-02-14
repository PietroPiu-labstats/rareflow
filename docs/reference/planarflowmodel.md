# Planar Normalizing Flow (1D)

A simple and readable implementation of a 1-dimensional planar flow: z_K
= z_0 + u \* h(wz_0 + b)

## Usage

``` r
planarflowmodel(u, w, b)
```

## Arguments

- u:

  Scalar parameter controlling the magnitude of the deformation.

- w:

  Scalar parameter controlling the slope of the activation.

- b:

  Scalar bias term.

## Value

A flow model object with methods:

- sampleq(n)

- logq(z0)

- applyflow(z0)

## Details

where:

- h is a smooth activation (tanh)

- the log-determinant is computed analytically

This flow is mainly useful for pedagogical purposes or as a lightweight
building block in variational inference.
