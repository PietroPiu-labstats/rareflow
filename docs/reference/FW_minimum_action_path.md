# Minimum Action Path (MAP)

Computes a simple Freidlin–Wentzell minimum action path between two
states using gradient descent on the discretized action functional.

## Usage

``` r
FW_minimum_action_path(
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

  Drift function b(x).

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

- path: matrix of size T × d

- action: FW action of the optimized path

## Details

This is a lightweight illustrative solver suitable for low-dimensional
systems. More advanced solvers (string method, MAM, etc.) can be
integrated.
