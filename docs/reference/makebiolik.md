# Biological Two-Separator Likelihood

Same structure as
[`makeneurolik()`](https://pietropiu-labstats.github.io/rareflow/reference/makeneurolik.md),
but with a different default separation parameter. This can be used to
model simple biological switching systems or coarse-grained gene
expression states.

## Usage

``` r
makebiolik(a = 0.2)
```

## Arguments

- a:

  Separation parameter controlling the spacing between the two logits.

## Value

A function mapping a latent vector `z` to a probability vector.
