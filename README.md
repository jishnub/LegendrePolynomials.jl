# LegendrePolynomials.jl

[![Build Status](https://travis-ci.com/jishnub/LegendrePolynomials.jl.svg?branch=master)](https://travis-ci.com/jishnub/LegendrePolynomials.jl)

Legendre polynomials and their derivatives using automatic differentiation

# Getting Started

## Prerequisites

The package depends on `HyperDualNumbers` and `OffsetArrays`. 

## Installing

To install the package, run 

```julia
] add LegendrePolynomials
```

## Using

This package computes Legendre polynomials and their first and second derivatives in one pass using hyperdual numbers. All the functions require an argument `x` that satisfies `-1<=x<=1`, and either an integer `l` that indicates the degree of the polynomial, or  a possibly optional keyword argument `lmax` that indicates the maximum degree of the Legendre polynomials that will be computed. The polynomials are computed iteratively using a three-term recurrence relation, so all the values in the range will have to be computed.

There are two classes of functions, one that returns the value of the Legendre polynomial and its derivatives for a single `l`, and one that returns the values for all `l` in `0:lmax` for a specified `lmax` as an array. The arrays returned will have the Legendre polynomials and their derivatives stored along columns.

There are three quantities that can be computed: `Pl`, `dPl` and `d2Pl`. It's possible to compute them individually or as combinations of two or three at once. There are allocating as well as non-allocating methods that do the same.

### Individual degrees

The way to compute Legendre polynomials of degree `l` and argument `x` is through the function `Pl(x,l)`. This is the general signature followed by the other functions as well that compute the derivatives for a single degree. The functions are relatively accurate and performant till degrees of millions, although a thorough test of accuracy has not been carried out (verified roughly with results from Mathematica for some degrees).

```julia
julia> Pl(0.5,10)
-0.18822860717773438

julia> @btime Pl(0.5,1_000_000)
  13.022 ms (0 allocations: 0 bytes)
-0.0006062610545162491
```

The accuracy can be increased by using Arbitrary Precision Arithmetic, through the use of 
`BigInt` and `BigFloat`. This comes at the expense of performance though.

```julia
julia> @btime Pl(1/3,1_000)
  11.135 μs (0 allocations: 0 bytes)
0.01961873093750127

julia> @btime Pl(big(1)/3,1_000)
  1.751 ms (58965 allocations: 3.13 MiB)
0.01961873093750094969323575593064773353450010511010742490834078609343359408691498
```

The way to compute the derivatives is through `dPl(x,l)` and `d2Pl(x,l)`, for example

```julia
julia> dPl(0.5,5)
-2.2265625000000004

julia> d2Pl(0.5,5)
-6.5625
```

Combinations of these three can be computed as well, for example

```julia
julia> Pl_dPl(0.5,5)
(0.08984375, -2.2265625000000004)

julia> Pl_d2Pl(0.5,5)
(0.08984375, -6.5625)

julia> dPl_d2Pl(0.5,5)
(-2.2265625000000004, -6.5625)

julia> Pl_dPl_d2Pl(0.5,5)
(0.08984375, -2.2265625000000004, -6.5625)
```

### All degrees up to a cutoff

The second class of methods return all the polynomials up to a cutoff degree `lmax`. They are returned as `OffsetArrays` that have 0-based indexing, keeping in mind that the polynomials start from `l=0`.

The polynomials and their derivatives can be computed in general by calling the function `P(x;lmax)`, where `P` has to be chosen appropriately as necessary. The keyword argument has to be specified in the allocating functions, whereas it may be inferred from the array in the non-allocating versions.

An example of the allocating functions are

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

It's possible to compute combinations analogous to the ones seen earlier, for example

```julia
julia> Pl_dPl_d2Pl(0.5,lmax=3)
([1.0, 0.5, -0.125, -0.4375], [0.0, 1.0, 1.5, 0.3750000000000001], [0.0, 0.0, 3.0, 7.5])
```

This returns a 3-tuple of `OffsetArrays` `Pl`, `dPl` and `d2Pl`. 

There are non-allocating functions as well that can be called using a pre-allocated array. As an example

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

julia> fill!(P,0);

julia> dPl!(P,0.5)
4×3 OffsetArray(::Array{Float64,2}, 0:3, 0:2) with eltype Float64 with indices 0:3×0:2:
 0.0  0.0    0.0
 0.0  1.0    0.0
 0.0  1.5    0.0
 0.0  0.375  0.0
```

Note that the column number that will be populated depends on the order of the derivative, assuming 0-based indexing. Therefore `Pl!` will fill column 0, `dPl!` will fill column 1 and `d2Pl!` will fill column 2. Combinations of these will fill multiple columns as expected.

# License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/jishnub/LegendrePolynomials.jl/blob/master/LICENSE) file for details.