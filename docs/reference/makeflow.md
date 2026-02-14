# Flow Factory

A unified constructor for all flow models in the rareflow package. This
function dispatches to the appropriate flow implementation based on the
`flowtype` argument.

## Usage

``` r
makeflow(flowtype, params = list())
```

## Arguments

- flowtype:

  Character string specifying the flow type.

- params:

  A named list of parameters required by the chosen flow.

## Value

A flow model object with methods:

- sampleq(n)

- logq(z0)

- applyflow(z0)

## Details

Supported flow types:

- "planar" -\> planarflowmodel()

- "radial" -\> radialflowmodel()

- "maf" -\> mafflowmodel()

- "splinepwlin" -\> splinepwlinflowmodel()

## Examples

``` r
# Create a planar flow
f <- makeflow("planar", list(u = 0.1, w = 0.2, b = 0))
s <- f$sampleq(10)

# Create a 2D spline flow
f2 <- makeflow("splinepwlin", list(d = 2, K = 8))
```
