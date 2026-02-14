# rareflow: Variational Flow-Based Inference for Rare Events

**Normalizing Flows for Rare-Event Inference**

![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)

License: GPL-3

------------------------------------------------------------------------

## Overview

**rareflow** provides a unified framework for rare-event inference by
combining:

- **Sanov theory** for empirical distributions  
- **Girsanov change of measure** for SDEs  
- **Freidlin–Wentzell large deviations** for small-noise diffusions  
- **Normalizing flows** for flexible variational inference

The package includes:

- modular flow models (planar, radial, MAF, spline)  
- variational optimization via ELBO  
- wrappers for Girsanov tilting and FW quasi-potentials  
- tools for rare-event simulation and analysis

------------------------------------------------------------------------

## Installation

### From GitHub

``` r
# install.packages("devtools")
devtools::install_github("PietroPiu-labstats/rareflow")
```

## Quick Example

A minimal workflow for fitting a variational posterior using a planar
flow:

``` r
library(rareflow)

Qobs <- c(0.05, 0.90, 0.05)
px <- function(z) c(0.3, 0.4, 0.3)

flow <- makeflow("planar")
fit <- fitflowvariational(Qobs, pxgivenz = px, nmc = 500)

fit$elbo
```

## Girsanov Tilting

``` r
b <- function(x) -x
dt <- 0.01
T <- 1000
Winc <- rnorm(T, sd = sqrt(dt))

fit_gir <- fitflow_girsanov(
  observed = Qobs,
  drift = b,
  Winc = Winc,
  dt = dt,
  pxgivenz = px
)
```

## Freidlin–Wentzell Quasi-Potential

``` r
b <- function(x) x - x^3
qp <- FW_quasipotential(-1, 1, drift = b, T = 200, dt = 0.01)
qp$action
```

## Documentation

Full documentation and examples are available in the package vignette:

``` r
vignette("rareflow")
```

## Features

- Variational inference with normalizing flows

- Rare-event tilting via Girsanov

- Minimum-action paths and quasi-potentials

- Support for 1D and 2D systems

- Modular flow architecture
