# Freidlin-Wentzell Quasi-Potential via Path Minimization

Computes an approximate Freidlin-Wentzell quasi-potential between two
points \\x_0\\ and \\x_1\\ by minimizing the FW action functional over
discretized paths.

## Usage

``` r
FW_quasipotential(
  x0,
  x1,
  drift,
  T = 200,
  dt = 0.01,
  niter = 200,
  stepsize = 0.1
)
```

## Arguments

- x0:

  Starting point (numeric vector).

- x1:

  Target point (numeric vector).

- drift:

  Drift function \\b(x)\\.

- T:

  Number of time steps.

- dt:

  Time step.

- niter:

  Number of gradient descent iterations.

- stepsize:

  Gradient descent step size.

## Value

A list with:

- path: matrix of size \\T \times d\\

- action: FW action of the optimized path

## Details

The algorithm:

1.  Initializes a straight-line path between \\x_0\\ and \\x_1\\.

2.  Performs simple gradient descent on the FW action.

This is a naive but effective illustrative method for low-dimensional
systems. More advanced solvers (string method, MAM, etc.) can be plugged
in.
