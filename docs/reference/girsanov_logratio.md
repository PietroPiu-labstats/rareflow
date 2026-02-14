# Girsanov Log-Ratio for Drift-Tilted Diffusions

Computes the Radon–Nikodym derivative (log form) associated with a
Girsanov change of measure for an SDE:

## Usage

``` r
girsanov_logratio(theta_path, Winc, dt)
```

## Arguments

- theta_path:

  Numeric vector of drift tilts \\\theta_t\\.

- Winc:

  Numeric vector of Brownian increments \\\Delta W_t\\.

- dt:

  Time step size.

## Value

A numeric log-likelihood ratio.

## Details

\\dX_t = b(X_t)\\ dt + dW_t\\

tilted by an alternative drift:

\\dX_t = (b(X_t) + \theta_t)\\ dt + dW_t\\.

The log-likelihood ratio is:

\\\log \frac{dQ}{dP} = \sum_t \left( \theta_t W\_{t} -
\frac{1}{2}\theta_t^2\\ dt \right)\\.

This function returns the log-ratio for a given path of drift tilts
`theta_path`, Brownian increments `Winc`, and time step `dt`.
