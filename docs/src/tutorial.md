```@meta
DocTestSetup = :(using LegendrePolynomials)
```

## Computing normalized Legendre polynomials

The norm of the polynomials may be specified through the keyword argument `norm`. The various normalization options are listed in the [Normalization options](@ref) section.

```jldoctest
julia> l = 3;

julia> Pl(0.5, l)
-0.4375

julia> Pl(0.5, l, norm = Val(:normalized))
-0.8184875533567997

julia> Pl(0.5, l, norm = Val(:normalized)) ≈ √((2l+1)/2) * Pl(0.5, l)
true

julia> l = m = 3000;

julia> Plm(0.5, l, m, norm = Val(:normalized))
2.172276347346834e-187
```

Starting from these, other normalization such as the [Ambix SN3D format](https://en.wikipedia.org/wiki/Ambisonic_data_exchange_formats#SN3D) may be constructed as

```julia
julia> Plmambix(x, l, m) = Plm(x, l, m, norm=Val(:schmidtquasi)) / √(4pi)
Plmambix (generic function with 1 method)

julia> Plmambix(0.5, 2, 1)
-0.21157109383040862
```

## Condon-Shortley phase

By default, the associated Legendre polynomials are computed by including the Condon-Shortley phase ``(-1)^m``.
This may be disabled by passing the flag `csphase = false` to the various `Plm` functions.

!!! note
    The `Pl` functions do not accept the keyword argument `csphase`.

```jldoctest
julia> Plm(0.5, 3, 3)
-9.742785792574933

julia> Plm(0.5, 3, 3, csphase = false)
9.742785792574933

julia> all(Plm(0.5, 3, m, csphase = false) ≈ (-1)^m * Plm(0.5, 3, m) for m in -3:3)
true
```

## Polynomials at multiple points

The `collectPl(m)` functions return a vector of polynomials evaluated for one ``\theta``. To evaluate the polynomials for multiple ``\theta`` in one go, one may run

```jldoctest multipletheta
julia> θr = range(0, pi, length=6);

julia> collectPl.(cos.(θr), lmax=3)
6-element Vector{OffsetArrays.OffsetVector{Float64, Vector{Float64}}}:
 [1.0, 1.0, 1.0, 1.0]
 [1.0, 0.8090169943749475, 0.4817627457812106, 0.1102457514062632]
 [1.0, 0.30901699437494745, -0.35676274578121053, -0.38975424859373686]
 [1.0, -0.30901699437494734, -0.35676274578121064, 0.3897542485937368]
 [1.0, -0.8090169943749473, 0.48176274578121037, -0.11024575140626285]
 [1.0, -1.0, 1.0, -1.0]
```

This returns a vector of vectors, where each element corresponds to one ``\theta``. Often, one wants a `Matrix` where
each column corresponds to one ``\theta``. We may obtain this as

```jldoctest multipletheta
julia> mapreduce(hcat, θr) do θ
           collectPl(cos(θ), lmax=3)
       end
4×6 Matrix{Float64}:
 1.0  1.0        1.0        1.0        1.0        1.0
 1.0  0.809017   0.309017  -0.309017  -0.809017  -1.0
 1.0  0.481763  -0.356763  -0.356763   0.481763   1.0
 1.0  0.110246  -0.389754   0.389754  -0.110246  -1.0
```

As of Julia `v1.7.2` and `OffsetArrays` `v1.10.8`, this strips off the axis information, so one would need to wrap the result in an `OffsetArray` again to index the `Matrix` using the degrees ``\ell``.

## Increasing precision

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

This is particularly important to avoid overflow while computing high-order derivatives. For example:

```jldoctest
julia> dnPl(0.5, 300, 200) # Float64
NaN

julia> dnPl(big(1)/2, 300, 200) # BigFloat
1.738632750542319394663553898425873258768856732308227932150592526951212145232716e+499
```

In general, one would need to use higher precision for both the argument `x` and the degrees `l` to obtain accurate results. For example:

```jldoctest
julia> Plm(big(1)/2, big(3000), 3000)
2.05451584347939475644802290993338963448971107938391335496027846832433343889916e+9844
```

## Symbolic evaluation

It's possible to symbolically evaluate Legendre or associated Legendre polynomials using [`Symbolics.jl`](https://github.com/JuliaSymbolics/Symbolics.jl). An exmaple is:

```jldoctest
julia> using Symbolics

julia> @variables x;

julia> Pl(x, 3)
(5//3)*x*((3//2)*(x^2) - (1//2)) - (2//3)*x

julia> myf = eval(build_function(Pl(x, 3), [x]));

julia> myf(0.4) == Pl(0.4, 3)
true
```
