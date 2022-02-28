module LegendrePolynomials

using OffsetArrays
using SpecialFunctions

export Pl
export collectPl
export collectPl!
export Plm
export collectPlm
export collectPlm!
export dnPl
export collectdnPl
export collectdnPl!

function checkdomain(x)
    abs(x) > 1 && throw(DomainError(x,
        "Legendre Polynomials are defined for arguments lying in -1 ⩽ x ⩽ 1"))
end

assertnonnegative(l) = l >= 0 || throw(ArgumentError("l must be >= 0, received " * string(l)))
assertnonnegative(::Nothing) = nothing

checklm(l::Integer, m::Integer) = 0 <= m <= l || throw(ArgumentError("m must satisfy 0 <= m <= l"))
checklm(::Any, ::Any) = nothing

function checklength(arr, minlength)
    length(arr) >= minlength || throw(ArgumentError(
        "array is not large enough to store all values, require a minimum length of " * string(minlength)))
end

@inline function Pl_recursion(::Type{T}, ℓ, (Plm1, Plm2, α), x) where {T}
  # relation is valid from ℓ = 1
  Pl = ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ
  convert(T, Pl), α
end

@inline function Plm_recursion(::Type{T}, ℓ, m, (Plm1, Plm2, α_ℓm1_m), x) where {T}
    # relation is valid from ℓ = m+2
    α_ℓ_m = coeff(ℓ, m)
    Pl = α_ℓ_m * (x * Plm1 - inv(α_ℓm1_m) * Plm2)
    return convert(T, Pl), α_ℓ_m
end

@inline function dPl_recursion(::Type{T}, ℓ, n, P_n_lm1, P_nm1_lm1, P_n_lm2, x) where {T}
  P_n_l = ((2ℓ - 1) * (x * P_n_lm1 + n * P_nm1_lm1) - (ℓ - 1) * P_n_lm2)/ℓ
  convert(T, P_n_l)
end

# special cases
# case 1 : l == n, in which case there are no (l-1,n) and (l-2,n) terms
@inline function dPl_recursion(::Type{T}, ℓ, n, ::Nothing, P_nm1_lm1, ::Nothing, x) where {T}
    P_n_l = (2ℓ-1) * P_nm1_lm1
    convert(T, P_n_l)
end
# case 1 : l == n + 1, in which case there's no (l-2,n) term
@inline function dPl_recursion(::Type{T}, ℓ, n, P_n_lm1, P_nm1_lm1, ::Nothing, x) where {T}
    P_n_l = ((2ℓ - 1) * (x * P_n_lm1 + n * P_nm1_lm1))/ℓ
    convert(T, P_n_l)
end

function logplm_norm(l, m)
    T = promote_type(typeof(l), typeof(m))
    (log(T(2)) - log(2T(l)+1) + logfactorial(l + m) - logfactorial(l - m))/2
end
function _maybebigexp(t)
    if t < log(prevfloat(typemax(t)))
        return exp(t)
    else
        return exp(big(t))
    end
end
function plm_norm(l, m)
    t = logplm_norm(l, m)
    _maybebigexp(t)
end
normtype(l, m, x, ::Val{:normalized}) = polytype(x)
normtype(l, m, x, ::Val{:standard}) = typeof(plm_norm(l, m))
function logabspll_prefactor(l, T = typeof(l))
    lT = T(l)
    a = (logfactorial(2lT-1) + log(2lT+1) - log(lT))/2
    b = l*log(T(2)) + logfactorial(lT-1)
    a - b
end
neg1pow(l) = iseven(l) ? 1 : -1
function pll_prefactor(l, T = typeof(l))
    l == 0 && return sqrt(oftype(logfactorial(T(l)), 1/2))
    t = logabspll_prefactor(l, T)
    neg1pow(l) * exp(t)
end
function coeff(l, m, T = float(promote_type(typeof(l), typeof(m))))
    √((2T(l) - 1)/((T(l)-m)*(T(l)+m))*(2T(l)+1))
end

polytype(x) = float(typeof(x))

struct LegendrePolynomialIterator{T, L<:Union{Integer,Nothing}, M<:Union{Integer, Nothing}, V}
    x :: V
    lmax :: L
    m :: M
    function LegendrePolynomialIterator{T,L,M,V}(x::V, lmax::L, m::M) where {T, L<:Union{Integer,Nothing}, M<:Union{Integer, Nothing}, V}
        checkdomain(x)
        new{T,L,M,V}(x, lmax, m)
    end
end
LegendrePolynomialIterator(x) = LegendrePolynomialIterator{polytype(x), Nothing, Nothing, typeof(x)}(x, nothing, nothing)
function LegendrePolynomialIterator(x, lmax::Integer)
    assertnonnegative(lmax)
    LegendrePolynomialIterator{polytype(x), typeof(lmax), Nothing, typeof(x)}(x, lmax, nothing)
end
function LegendrePolynomialIterator(x, lmax::Integer, m::Integer)
    assertnonnegative(lmax)
    checklm(lmax, m)
    lmaxT, mT = promote(lmax, m)
    LegendrePolynomialIterator{polytype(x), typeof(lmaxT), typeof(mT), typeof(x)}(x, lmaxT, mT)
end

Base.eltype(::Type{<:LegendrePolynomialIterator{T}}) where {T} = T
Base.IteratorSize(::Type{<:LegendrePolynomialIterator{<:Any,Nothing}}) = Base.IsInfinite()

function Base.iterate(iter::LegendrePolynomialIterator{T,<:Any,<:Nothing}) where {T}
    Pl = one(T)
    Plm1 = zero(T)
    nextstate = (Pl, Plm1, nothing)
    return Pl, (1, nextstate)
end
# starting value, l == m
function Base.iterate(iter::LegendrePolynomialIterator{T,<:Any,<:Integer}) where {T}
    m = iter.m
    if m == 0
        Pl = one(T)
    else
        x = iter.x
        f = pll_prefactor(m)
        t = f * (√(1-x^2))^m
        Pl = convert(T, t)
    end
    Plm1 = zero(T)
    α_ℓ_m = Inf
    return Pl, (m+1, (Pl, Plm1, α_ℓ_m))
end
_isdone(l, lmax) = l > lmax
_isdone(::Any, ::Nothing) = false

function Base.iterate(iter::LegendrePolynomialIterator{T,<:Any,Nothing}, state) where {T}
    l, prevstate = state
    _isdone(l, iter.lmax) && return nothing
    x = iter.x
    # standard-norm Bonnet iteration
    Pl, α = Pl_recursion(T, l, prevstate, x)
    Plm1 = first(prevstate)
    nextstate = (Pl, Plm1, α)
    return Pl, (l+1, nextstate)
end
function Base.iterate(iter::LegendrePolynomialIterator{T,<:Any,<:Integer}, state) where {T}
    l, prevstate = state
    _isdone(l, iter.lmax) && return nothing
    m = iter.m
    x = iter.x
    if m == 0
        # standard-norm Bonnet iteration
        Pl, α_ℓ_m = Pl_recursion(T, l, prevstate, x)
    else
        # iterate over normalized functions
        Pl, α_ℓ_m = Plm_recursion(T, l, m, prevstate, x)
    end
    Plm1 = first(prevstate)
    nextstate = (Pl, Plm1, α_ℓ_m)
    return Pl, (l+1, nextstate)
end

_length(iter) = _length(iter.lmax, iter.m)
_length(lmax::Integer, m::Integer) = lmax - m + 1
_length(lmax::Integer, m::Nothing) = lmax + 1
_axes1(iter) = _axes1(iter.lmax, iter.m)
_axes1(lmax::Integer, m::Integer) = m:lmax
_axes1(lmax::Integer, m::Nothing) = 0:lmax
Base.length(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = _length(iter)
Base.size(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = (length(iter),)
Base.axes(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = (_axes1(iter), )
Base.keys(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = _axes1(iter)

Base.copy(iter::LegendrePolynomialIterator) = typeof(iter)(iter.x, iter.lmax, iter.m)

maybenormalize(P, l, ::Val{:normalized}) = oftype(P, P * √(2/(2l+1)))
maybenormalize(P, l, ::Val{:standard}) = P

function maybenormalize(P, l, m, norm::Val{:normalized})
    if m == 0
        return maybenormalize(P, l, norm)
    else
        return P
    end
end
function maybenormalize(P, l, m, norm::Val{:standard})
    if m == 0
        return maybenormalize(P, l, norm)
    else
        return plm_norm(l, m) * P
    end
end
maybenormalize(P, l, m, norm::Val) = throw(ArgumentError("norm = $norm undefined, valid norms are :standard and :normalized"))

"""
    Pl(x, l::Integer; [norm = Val(:standard)])

Compute the Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and the degree `l`.

The default norm is chosen to be `Val(:standard)`, in which case the polynomials
satisfy `Pl(1, l) == 1` for all `l`.
Optionally, the `norm` may be set to `Val(:normalized)` to obtain normalized Legendre polynomials.
These have an L2 norm of `1`.

# Examples
```jldoctest
julia> Pl(1, 2)
1.0

julia> Pl(0.5, 4) ≈ -37/128 # analytically obtained value
true

julia> Pl(0.5, 20, norm = Val(:normalized))
-0.010680579639558057
```
"""
function Pl(x, l::Integer; norm = Val(:standard))
    _checkvalues(x, l, 0, 0)
    T = polytype(x)
    if x == 1
        return one(T)
    elseif x == -1
        return convert(T, neg1pow(l))
    elseif x == 0 && isodd(l)
        return zero(T)
    end
    iter_full = LegendrePolynomialIterator(x, l)
    iter = Iterators.drop(iter_full, length(iter_full) - 1)
    Pl = only(iter)
    return maybenormalize(Pl, l, norm)
end

function doublefactorial(T, n)
    p = convert(T, one(n))
    for i in n:-2:1
        p *= convert(T, i)
    end
    convert(T, p)
end

"""
    Plm(x, l::Integer, m::Integer; [norm = Val(:standard)])

Compute the associated Legendre polynomial ``P_\\ell^m(x)`` for a non-negative ``m``.
Optionally, the `norm` may be set to `Val(:normalized)`, in which case normalized
associated Legendre polynomials are evaluated. These have an L2 norm of `1`.

The coefficient `m` must be non-negative. For `m == 0` this function just returns
Legendre polynomials.

# Examples

```jldoctest
julia> Plm(0.5, 3, 2) ≈ 45/8 # analytically obtained value
true

julia> Plm(0.5, 4, 0) == Pl(0.5, 4)
true
```
"""
function Plm(x, l::Integer, m::Integer; norm = Val(:standard))
    _checkvalues(x, l, m, 0)
    if l < m
        return zero(polytype(x))
    end
    # special cases
    T = polytype(x)

    lT, mT = promote(l, m)

    iter_full = LegendrePolynomialIterator(x, lT, mT)
    iter = Iterators.drop(iter_full, length(iter_full) - 1)
    Plm = only(iter)
    return maybenormalize(Plm, lT, mT, norm)
end

function _checkvalues(x, l, m, n)
    checkdomain(x)
    assertnonnegative(l)
    checklm(l, m)
    n >= 0 || throw(ArgumentError("order of derivative n must be >= 0"))
end

Base.@propagate_inbounds function _unsafednPl!(cache, x, l, n)
    # unsafe, assumes 1-based indexing
    checklength(cache, l - n + 1)
    if n == l # may short-circuit this
        cache[1] = doublefactorial(eltype(cache), 2l-1)
    else
        collectPl!(cache, x, lmax = l - n)

        for ni in 1:n
            # We denote the terms as P_ni_li

            # li == ni
            P_nim1_nim1 = cache[1]
            P_ni_ni = dPl_recursion(eltype(cache), ni, ni, nothing, P_nim1_nim1, nothing, x)
            cache[1] = P_ni_ni

            # li == ni + 1
            P_nim1_ni = cache[2]
            P_ni_nip1 = dPl_recursion(eltype(cache), ni + 1, ni, P_ni_ni, P_nim1_ni, nothing, x)
            cache[2] = P_ni_nip1

            for li in ni+2:min(l, l - n + ni)
                P_ni_lim2 = cache[li - ni - 1]
                P_ni_lim1 = cache[li - ni]
                P_nim1_lim1 = cache[li - ni + 1]
                P_ni_li = dPl_recursion(eltype(cache), li, ni, P_ni_lim1, P_nim1_lim1, P_ni_lim2, x)
                cache[li - ni + 1] = P_ni_li
            end
        end
    end
    nothing
end

"""
    dnPl(x, l::Integer, n::Integer, [cache::AbstractVector])

Compute the ``n``-th derivative ``d^n P_\\ell(x)/dx^n`` of the Legendre polynomial ``P_\\ell(x)``.
Optionally a pre-allocated vector `cache` may be provided, which must have a minimum length of `l - n + 1`
and may be overwritten during the computation.

The order of the derivative `n` must be non-negative. For `n == 0` this function just returns
Legendre polynomials.

# Examples

```jldoctest
julia> dnPl(0.5, 3, 2) # second derivative of P3(x) at x = 0.5
7.5

julia> dnPl(0.5, 4, 0) == Pl(0.5, 4) # zero-th order derivative == Pl(x)
true
```
"""
Base.@propagate_inbounds function dnPl(x, l::Integer, n::Integer,
    A = begin
        _checkvalues(x, l, 0, n)
        # do not allocate A if the value is trivially zero
        if l < n
            return zero(polytype(x))
        end
        zeros(polytype(x), l - n + 1)
    end
    )

    _checkvalues(x, l, 0, n)
    # check if the value is trivially zero in case A is provided in the function call
    if l < n
        return zero(eltype(A))
    end

    cache = OffsetArrays.no_offset_view(A)
    # function barrier, as no_offset_view may be type-unstable
    _unsafednPl!(cache, x, l, n)

    return cache[l - n + 1]
end

"""
    collectPl!(v::AbstractVector, x; [lmax::Integer = length(v) - 1], [norm = Val(:standard)])

Compute the Legendre Polynomials ``P_\\ell(x)`` for the argument `x` and all degrees `l` in `0:lmax`,
and store them in `v`.

At output `v[firstindex(v) + l] == Pl(x,l)`.

# Examples
```jldoctest
julia> v = zeros(4);

julia> collectPl!(v, 0.5);

julia> all(zip(0:3, v)) do (l, vl); vl ≈ Pl(0.5, l); end
true

julia> collectPl!(v, 0.5, norm = Val(:normalized));

julia> all(zip(0:3, v)) do (l,vl); vl ≈ Pl(0.5, l, norm = Val(:normalized)); end
true

julia> v = zeros(0:4);

julia> collectPl!(v, 0.5, lmax = 3) # only l from 0 to 3 are populated
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
  1.0
  0.5
 -0.125
 -0.4375
  0.0
```
"""
function collectPl!(v::AbstractVector, x; lmax::Integer  = length(v) - 1, norm = Val(:standard))
    _checkvalues(x, lmax, 0, 0)
    n = lmax + 1
    checklength(v, n)

    iter = LegendrePolynomialIterator(x, lmax)

    inds = (firstindex(v)-1) .+ (1:n)
    v_section = @view v[inds]

    for (ind, (l, Pl)) in enumerate(pairs(iter))
        v_section[ind] = maybenormalize(Pl, l, norm)
    end

    return v
end

"""
    collectPl(x; lmax::Integer, [norm = Val(:standard)])

Compute the Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and all degrees `l` in `0:lmax`.
Return `P` with indices `0:lmax`, where `P[l] == Pl(x,l)`

# Examples
```jldoctest
julia> v = collectPl(0.5, lmax = 4);

julia> all(l -> v[l] ≈ Pl(0.5, l), 0:4)
true

julia> v = collectPl(0.5, lmax = 4, norm = Val(:normalized));

julia> all(l -> v[l] ≈ Pl(0.5, l, norm = Val(:normalized)), 0:4)
true
```
"""
function collectPl(x; lmax::Integer, norm = Val(:standard))
    _checkvalues(x, lmax, 0, 0)
    v_offset = zeros(polytype(x), 0:lmax)
    v = parent(v_offset)
    collectPl!(v, x; lmax = lmax, norm = norm)
    return v_offset
end

"""
    collectPlm(x; lmax::Integer, m::Integer, [norm = Val(:standard)])

Compute the associated Legendre Polynomial ``P_\\ell,m(x)`` for the argument `x` and all degrees `l = 0:lmax`.

The polynomials are defined to include the Condon-Shortley phase ``(-1)^m``.

The coefficient `m` must be greater than or equal to zero.

Returns `v` with indices `m:lmax`, where `v[l] == Plm(x, l, m)`.

# Examples

```jldoctest
julia> v = collectPlm(0.5, lmax = 3, m = 2);

julia> all(l -> v[l] ≈ Plm(0.5, l, 2), 2:3)
true

julia> v = collectPlm(0.5, lmax = 3, m = 2, norm = Val(:normalized));

julia> all(l -> v[l] ≈ Plm(0.5, l, 2, norm = Val(:normalized)), 2:3)
true
```
"""
function collectPlm(x; lmax::Integer, m::Integer, norm = Val(:standard),
        Tnorm = begin
            norm === Val(:standard) ? begin
            _checkvalues(x, lmax, m, 0)
            typeof(plm_norm(lmax, m))
            end : float(promote_type(typeof(m), typeof(lmax)))
        end
        )
    _checkvalues(x, lmax, m, 0)
    T = promote_type(polytype(x), Tnorm)
    v = zeros(T, m:lmax)
    collectPlm!(parent(v), x; lmax = lmax, m = m, norm = norm)
    v
end

"""
    collectPlm!(v::AbstractVector, x; lmax::Integer, m::Integer, [norm = Val(:standard)])

Compute the associated Legendre Polynomial ``P_\\ell,m(x)`` for the argument `x` and all degrees `l = 0:lmax`,
and store the result in `v`.

The coefficient `m` must be greater than or equal to zero.

At output, `v[l + firstindex(v)] == Plm(x, l, m)` for `l = 0:lmax`.

# Examples

```jldoctest
julia> v = zeros(2);

julia> collectPlm!(v, 0.5, lmax = 3, m = 2);

julia> all(zip(2:3, v)) do (l, vl); vl ≈ Plm(0.5, l, 2); end
true
```
"""
function collectPlm!(v, x; m::Integer, lmax::Integer = length(v) - m + 1, norm = Val(:standard))
    assertnonnegative(lmax)
    m >= 0 || throw(ArgumentError("coefficient m must be >= 0"))
    checklm(lmax, m)
    n = lmax - m + 1
    checklength(v, n)

    iter = LegendrePolynomialIterator(x, lmax, m)

    inds = (firstindex(v)-1) .+ (1:n)
    v_section = @view v[inds]

    for (ind, (l, Plm)) in enumerate(pairs(iter))
        v_section[ind] = maybenormalize(Plm, l, m, norm)
    end

    return v
end

"""
    collectdnPl(x; lmax::Integer, n::Integer)

Compute the ``n``-th derivative of a Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and all degrees `l = 0:lmax`.

The order of the derivative `n` must be greater than or equal to zero.

Returns `v` with indices `0:lmax`, where `v[l] == dnPl(x, l, n)`.

# Examples

```jldoctest
julia> collectdnPl(0.5, lmax = 3, n = 2)
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
 0.0
 0.0
 3.0
 7.5
```
"""
function collectdnPl(x; lmax::Integer, n::Integer)
    assertnonnegative(lmax)
    n >= 0 || throw(ArgumentError("order of derivative n must be >= 0"))
    v = zeros(polytype(x), 0:lmax)
    if lmax >= n
        collectdnPl!(parent(v), x; lmax = lmax, n = n)
    end
    v
end

"""
    collectdnPl!(v::AbstractVector, x; lmax::Integer, n::Integer)

Compute the ``n``-th derivative of a Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and all degrees `l = 0:lmax`,
and store the result in `v`.

The order of the derivative `n` must be greater than or equal to zero.

At output, `v[l + firstindex(v)] == dnPl(x, l, n)` for `l = 0:lmax`.

# Examples

```jldoctest
julia> v = zeros(4);

julia> collectdnPl!(v, 0.5, lmax = 3, n = 2)
4-element Vector{Float64}:
 0.0
 0.0
 3.0
 7.5
```
"""
function collectdnPl!(v, x; lmax::Integer, n::Integer)
    assertnonnegative(lmax)
    n >= 0 || throw(ArgumentError("order of derivative n must be >= 0"))
    checklength(v, lmax + 1)

    # trivially zero for l < n
    fill!((@view v[(0:n-1) .+ firstindex(v)]), zero(eltype(v)))
    # populate the other elements
    @inbounds dnPl(x, lmax, n, @view v[(n:lmax) .+ firstindex(v)])

    v
end

end
