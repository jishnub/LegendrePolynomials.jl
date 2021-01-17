# LegendrePolynomials.jl

![CI](https://github.com/jishnub/LegendrePolynomials.jl/workflows/CI/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/jishnub/LegendrePolynomials.jl/branch/master/graph/badge.svg?token=PGMRgetYmz)](https://codecov.io/gh/jishnub/LegendrePolynomials.jl)
[![stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jishnub.github.io/LegendrePolynomials.jl/stable)
[![dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://jishnub.github.io/LegendrePolynomials.jl/dev)

Iterative computation of Legendre Polynomials

# Getting Started

## Installing

To install the package, run 

```julia
] add LegendrePolynomials
```

## Quick start

To compute the Legendre polynomials for a given argument `x` and a degree `l`, use `Pl(x,l)`:

```julia
julia> Pl(0.5, 20)
-0.04835838106737356
```

To compute the n-th derivative of the Legendre polynomial of degree `l` at the argument `x`, use `dnPl(x, l, n)`:

```julia
julia> dnPl(0.5, 10, 5)
-61760.91796875
```

To compute all the polynomials for `0 <= l <= lmax`, use `collectPl(x; lmax)`

```julia
julia> collectPl(0.5, lmax = 5)
6-element OffsetArray(::Array{Float64,1}, 0:5) with eltype Float64 with indices 0:5:
  1.0
  0.5
 -0.125
 -0.4375
 -0.2890625
  0.08984375
```

To compute all the n-th derivatives for `0 <= l <= lmax`, use `collectdnPl(x; n, lmax)`

```julia
julia> collectdnPl(0.5, lmax = 5, n = 3)
6-element OffsetArray(::Array{Float64,1}, 0:5) with eltype Float64 with indices 0:5:
  0.0
  0.0
  0.0
 15.0
 52.5
 65.625
```

Check the documentation for other usage.

# License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/jishnub/LegendrePolynomials.jl/blob/master/LICENSE) file for details.