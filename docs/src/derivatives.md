# Derivatives of Legendre Polynomials

The Julia automatic differentiation framework may be used to compute the derivatives of Legendre polynomials alongside their values. Since the defintions of the polynomials are completely general, they may be called with dual or hyperdual numbers as arguments to evaluate derivarives in one go. 
We demonstrate one example of this using the package [`HyperDualNumbers.jl`](https://github.com/JuliaDiff/HyperDualNumbers.jl) v4:

```@meta
DocTestSetup = quote
	using LegendrePolynomials
end
```

```jldoctest hyperdual
julia> x = 0.5;

julia> Pl(x, 3)
-0.4375

julia> using HyperDualNumbers

julia> xh = Hyper(x, one(x), one(x), zero(x));

julia> p = Pl(xh, 3)
-0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂
```

The Legendre polynomial ``P_\ell(x)`` may be obtained using 

```jldoctest hyperdual
julia> realpart(p)
-0.4375
```

The first derivative ``dP_\ell(x)/dx`` may be obtained as 

```jldoctest hyperdual
julia> ε₁part(p)
0.375
```

The second derivative ``d^2P_\ell(x)/dx^2`` may be obtained using 

```jldoctest hyperdual
julia> ε₁ε₂part(p)
7.5
```

Something similar may also be evaluated for all `l` iteratively. For example, the function `collectPl` evaluates Legendre polynomials for the `HyperDualNumber` argument for a range of `l`:

```jldoctest hyperdual
julia> collectPl(xh, lmax=4)
5-element OffsetArray(::Array{Hyper{Float64},1}, 0:4) with eltype Hyper{Float64} with indices 0:4:
                1.0 + 0.0ε₁ + 0.0ε₂ + 0.0ε₁ε₂
                0.5 + 1.0ε₁ + 1.0ε₂ + 0.0ε₁ε₂
             -0.125 + 1.5ε₁ + 1.5ε₂ + 3.0ε₁ε₂
        -0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂
 -0.2890625 - 1.5625ε₁ - 1.5625ε₂ + 5.625ε₁ε₂
```

We may extract the first derivatives by broadcasting the function `ε₁part` on the array as:

```jldoctest hyperdual
julia> ε₁part.(collectPl(xh, lmax=4))
5-element OffsetArray(::Array{Float64,1}, 0:4) with eltype Float64 with indices 0:4:
  0.0
  1.0
  1.5
  0.375
 -1.5625
```

Similarly the function `ε₁ε₂part` may be used to obtain the second derivatives. 

Several convenience functions to compute the derivatives of Legendre polynomials were available in `LegendrePolynomials` v0.2, but have been removed in v0.3. The users are encouraged to implement convenience functions to extract the derivatives as necessary. As an exmaple, we may compute the polynomials and their first and second derivatives together as

```jldoctest
julia> using HyperDualNumbers

julia> function Pl_dPl_d2Pl(x; lmax)
           xh = Hyper(x, one(x), one(x), zero(x))
           p = collectPl(xh, lmax = lmax)
           realpart.(p), ε₁part.(p), ε₁ε₂part.(p)
       end
Pl_dPl_d2Pl (generic function with 1 method)

julia> Pl_dPl_d2Pl(0.5, lmax = 3)
([1.0, 0.5, -0.125, -0.4375], [0.0, 1.0, 1.5, 0.375], [0.0, 0.0, 3.0, 7.5])
```
