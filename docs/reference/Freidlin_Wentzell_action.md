# Freidlin–Wentzell Action Functional

Computes the discrete Freidlin–Wentzell action for a path \\\phi(t)\\
represented as a matrix of size T × d. The continuous action is:

## Usage

``` r
Freidlin_Wentzell_action(phi, drift, dt)
```

## Arguments

- phi:

  Matrix of path values (T × d).

- drift:

  Drift function \\b(x)\\ returning a numeric vector.

- dt:

  Time step.

## Value

Numeric action value.

## Details

\$\$ I\[\phi\] = \frac{1}{2} \int_0^T \\ \dot{\phi}(t) - b(\phi(t)) \\^2
dt, \$\$

and the discrete approximation is:

\$\$ I \approx \frac{1}{2} \sum\_{t=1}^{T-1} \\ (\phi\_{t+1} -
\phi_t)/dt - b(\phi_t) \\^2 \\ dt. \$\$
