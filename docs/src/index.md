```@meta
CurrentModule = LegendrePolynomials
DocTestSetup = quote
	using LegendrePolynomials
end
```

# Introduction

Compute [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials) using a 3-term recursion relation (Bonnet’s recursion formula).

``P_\ell(x) = \left((2\ell-1) x P_{\ell-1}(x) - (\ell-1)P_{\ell - 2}(x)\right)/\ell``

Currently this package evaluates the standard polynomials that satisfy ``P_\ell(1) = 1``. These are normalized as 

```math
\int P_m(x) P_n(x) dx = \frac{2}{2n+1} \delta_{mn}.
```

There are two main functions: 

* [`Pl(x,l)`](@ref Pl): this evaluates the Legendre polynomial for a given degree `l` at the argument `x`. The argument needs to satisfy `-1 <= x <= 1`.
* [`collectPl(x; lmax)`](@ref collectPl): this evaluates all the polynomials for `l` lying in `0:lmax` at the argument `x`. As before the argument needs to lie in the domain of validity. Functionally this is equivalent to `Pl.(x, 0:lmax)`, except `collectPl` evaluates the result in one pass, and is therefore faster. There is also the in-place version [`collectPl!`](@ref) that uses a pre-allocated array.

# Quick Start

Evaluate the Legendre polynomial for one `l` at an argument`x` as `Pl(x, l)`:

```jldoctest
julia> Pl(0.5, 3)
-0.4375
```

Evaluate all the polynomials for `l` in `0:lmax` as 

```jldoctest
julia> collectPl(0.5, lmax = 3)
4-element OffsetArray(::Array{Float64,1}, 0:3) with eltype Float64 with indices 0:3:
  1.0
  0.5
 -0.125
 -0.4375
```

# Increase precision

The precision of the result may be changed by using arbitrary-precision types such as `BigFloat`. For example, using `Float64` arguments we obtain

```jldoctest
julia> Pl(1/3, 5)
0.33333333333333337
```

whereas using `BigFloat`, we obtain

```jldoctest
julia> Pl(big(1)/3, 5)
0.3333333333333333333333333333333333333333333333333333333333333333333333333333305
```

The precision of the latter may be altered using `setprecision`, as 

```jldoctest
julia> setprecision(300) do 
       Pl(big(1)/3, 5)
       end
0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333317
```

# Reference

```@autodocs
Modules = [LegendrePolynomials]
```
