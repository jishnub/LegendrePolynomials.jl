# LegendrePolynomials.jl

![CI](https://github.com/jishnub/LegendrePolynomials.jl/workflows/CI/badge.svg?branch=master)
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

To compute all the polynomials for `0 <= l <= lmax`, use 

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

Check the documentation for other usage.

# License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/jishnub/LegendrePolynomials.jl/blob/master/LICENSE) file for details.