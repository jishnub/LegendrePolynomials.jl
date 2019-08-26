# Legendre.jl

[![Build Status](https://travis-ci.com/jishnub/Legendre.jl.svg?branch=master)](https://travis-ci.com/jishnub/Legendre.jl)

Legendre polynomials and their derivatives using automatic differentiation

# Getting Started

## Prerequisites

The package depends on `HyperDualNumbers` and `OffsetArrays`. 

## Installing

To install the package, run 

```julia
] add https://github.com/jishnub/Legendre.jl.git
using Legendre
```

## Using

This package computes Legendre polynomials and their first and second derivatives in one pass using hyperdual numbers. All the functions require an argument `x` that satisfies `-1<=x<=1`, and a possibly optional keyword argument `lmax` that indicates the maximum degree of the Legendre polynomials that will be computed. The polynomials are computed iteratively, so all the values in `0:lmax` will have to be computed.

The arrays returned will have the polynomials stored along the first dimension, and possibly the various derivatives along the second dimension in the non-allocating versions. 

The way to obtain all these is 

```julia
julia> Pl_dPl_d2Pl(0.5,lmax=3)
([1.0, 0.5, -0.125, -0.4375], [0.0, 1.0, 1.5, 0.3750000000000001], [0.0, 0.0, 3.0, 7.5])
```

This returns a 3-tuple of `OffsetArrays` Pl, dPl and d2Pl. We can obtain these individually using 

```julia
julia> Pl(0.5,lmax=3)
4-element OffsetArray(::Array{Float64,1}, 0:3) with eltype Float64 with indices 0:3:
  1.0   
  0.5   
 -0.125 
 -0.4375

julia> dPl(0.5,lmax=3)
4-element OffsetArray(::Array{Float64,1}, 0:3) with eltype Float64 with indices 0:3:
 0.0               
 1.0               
 1.5               
 0.3750000000000001

julia> d2Pl(0.5,lmax=3)
4-element OffsetArray(::Array{Float64,1}, 0:3) with eltype Float64 with indices 0:3:
 0.0
 0.0
 3.0
 7.5
```

There are non-allocating functions as well that can be called as 
```julia
julia> P=zeros(0:3,0:2)
4×3 OffsetArray(::Array{Float64,2}, 0:3, 0:2) with eltype Float64 with indices 0:3×0:2:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

julia> Pl_dPl_d2Pl!(P,0.5)
4×3 OffsetArray(::Array{Float64,2}, 0:3, 0:2) with eltype Float64 with indices 0:3×0:2:
  1.0     0.0    0.0
  0.5     1.0    0.0
 -0.125   1.5    3.0
 -0.4375  0.375  7.5
```

# License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/jishnub/Legendre.jl/blob/master/LICENSE) file for details.