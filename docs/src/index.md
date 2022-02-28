```@meta
CurrentModule = LegendrePolynomials
DocTestSetup = quote
	using LegendrePolynomials
end
```

# Introduction

Compute [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials), [Associated Legendre polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials), and their derivatives using a 3-term recursion relation (Bonnet’s recursion formula).

```math
P_\ell(x) = \left((2\ell-1) x P_{\ell-1}(x) - (\ell-1)P_{\ell - 2}(x)\right)/\ell
```

By default this package evaluates the standard polynomials that satisfy ``P_\ell(1) = 1`` and ``P_0(x) = 1``. These are normalized as

```math
\int_{-1}^1 P_m(x) P_n(x) dx = \frac{2}{2n+1} \delta_{mn}.
```

Optionally, normalized polynomials may be evaluated that have an L2 norm of `1`.

Analogous to Legendre polynomials, one may evaluate associated Legendre polynomials using a 3-term recursion relation. This is evaluated by iterating over the normalized associated Legendre functions, and multiplying the norm at the final stage. Such an iteration avoids floating-point overflow.

The relation used in evaluating associated Legendre polynomials is

```math
\bar{P}_{\ell}^{m}\left(x\right)=\alpha_{\ell m}\left(x\bar{P}_{\ell-1}^{m}\left(x\right)-\frac{1}{\alpha_{\ell-1m}}\bar{P}_{\ell-2}^{m}\left(x\right)\right),
```

where
```math
\alpha_{\ell m}=\sqrt{\frac{\left(2\ell+1\right)\left(2\ell-1\right)}{\left(\ell-m\right)\left(\ell+m\right)}},
```

There are six main functions:

* [`Pl(x,l)`](@ref Pl): this evaluates the Legendre polynomial for a given degree `l` at the argument `x`. The argument needs to satisfy `-1 <= x <= 1`.
* [`collectPl(x; lmax)`](@ref collectPl): this evaluates all the polynomials for `l` lying in `0:lmax` at the argument `x`. As before the argument needs to lie in the domain of validity. Functionally this is equivalent to `Pl.(x, 0:lmax)`, except `collectPl` evaluates the result in one pass, and is therefore faster. There is also the in-place version [`collectPl!`](@ref) that uses a pre-allocated array.
* [`Plm(x, l, m)`](@ref Plm): this evaluates the associated Legendre polynomial ``P_\ell,m(x)`` at the argument ``x``. The argument needs to satisfy `-1 <= x <= 1`.
* [`collectPlm(x; m, lmax)`](@ref collectPlm): this evaluates the associated Legendre polynomials with coefficient `m` for `l = 0:lmax`. There is also an in-place version [`collectPlm!`](@ref) that uses a pre-allocated array.
* [`dnPl(x, l, n)`](@ref dnPl): this evaluates the ``n``-th derivative of the Legendre polynomial ``P_\ell(x)`` at the argument ``x``. The argument needs to satisfy `-1 <= x <= 1`.
* [`collectdnPl(x; n, lmax)`](@ref collectdnPl): this evaluates the ``n``-th derivative of all the Legendre polynomials for `l = 0:lmax`. There is also an in-place version [`collectdnPl!`](@ref) that uses a pre-allocated array.

# Quick Start

Evaluate the Legendre polynomial for one `l` at an argument`x` as `Pl(x, l)`:

```jldoctest
julia> p = Pl(0.5, 3)
-0.4375

julia> p ≈ -7/16 # analytical value
true
```

Evaluate the associated Legendre Polynomial one `l,m` pair as `Plm(x, l, m)`:

```jldoctest
julia> p = Plm(0.5, 3, 2)
5.624999999999997

julia> p ≈ 45/8 # analytical value
true
```

Evaluate the `n`th derivative for one `l` as `dnPl(x, l, n)`:

```jldoctest
julia> dnPl(0.5, 3, 2)
7.5
```

Evaluate all the polynomials for `l` in `0:lmax` as `collectPl(x; lmax)`

```jldoctest
julia> collectPl(0.5, lmax = 3)
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
  1.0
  0.5
 -0.125
 -0.4375
```

Evaluate all the associated Legendre Polynomials for coefficient `m` as `collectPlm(x; lmax, m)`:

```jldoctest
julia> collectPlm(0.5, lmax = 5, m = 3)
3-element OffsetArray(::Vector{Float64}, 3:5) with eltype Float64 with indices 3:5:
  -9.742785792574933
 -34.09975027401223
 -42.62468784251535
```

Evaluate all the `n`th derivatives as `collectdnPl(x; lmax, n)`:

```jldoctest
julia> collectdnPl(0.5, lmax = 5, n = 3)
6-element OffsetArray(::Vector{Float64}, 0:5) with eltype Float64 with indices 0:5:
  0.0
  0.0
  0.0
 15.0
 52.5
 65.625
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

This is particularly important to avoid overflow while computing high-order derivatives. For example:

```jldoctest
julia> dnPl(0.5, 300, 200) # Float64
NaN

julia> dnPl(big(1)/2, 300, 200) # BigFloat
1.738632750542319394663553898425873258768856732308227932150592526951212145232716e+499
```

# Reference

```@autodocs
Modules = [LegendrePolynomials]
```
