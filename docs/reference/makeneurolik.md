# Neural Two-Separator Likelihood

Constructs a simple 3-category likelihood model based on a latent vector
`z`. The likelihood is defined by two logistic separators:

## Usage

``` r
makeneurolik(a = 0.3)
```

## Arguments

- a:

  Separation parameter controlling the spacing between the two logits.

## Value

A function mapping a latent vector `z` to a probability vector.

## Details

p1 = sigmoid(mean(z) - a) p2 = sigmoid(mean(z) + a)

producing a 3-class probability vector:

(1 - p1, p1 - p2, p2)

This likelihood is useful for toy neural classification models or simple
latent-to-categorical mappings.
