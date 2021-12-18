module LegendrePolynomials

using OffsetArrays

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

assertnonnegative(l) = (l >= 0 || throw(ArgumentError("l must be >= 0, received " * string(l))))

function checklength(arr, minlength)
	length(arr) >= minlength || throw(ArgumentError(
		"array is not large enough to store all values, require a minimum length of " * string(minlength)))
end

@inline function Pl_recursion(::Type{T}, ℓ, Plm1, Plm2, x) where {T}
	# relation is valid from ℓ = 1
	Pl = ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ
	convert(T, Pl)
end

@inline function Plm_recursion(::Type{T}, ℓ, m, P_m_lm1, P_m_lm2, x) where {T}
	P_m_l = ((2ℓ-1) * x * P_m_lm1 - (ℓ+m-1) * P_m_lm2) / (ℓ-m)
	convert(T, P_m_l)
end

# special case
# l == m, in which case there are no (l-1,m) and (l-2,m) terms
@inline function Plm_recursion_m(::Type{T}, ℓ, m, P_mm1_lm1, x) where {T}
	P_m_l = -(2ℓ-1) * sqrt(1 - x^2) * P_mm1_lm1
	convert(T, P_m_l)
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

polytype(x) = typeof(float(x*x))

"""
	LegendrePolynomialIterator(x, [lmax::Integer])

Return an iterator that generates the values of the Legendre polynomials ``P_\\ell(x)`` for the given `x`.
If `lmax` is specified then only the values of ``P_\\ell(x)`` from `0` to `lmax` are returned.

# Examples
``jldoctest
julia> import LegendrePolynomials: LegendrePolynomialIterator

julia> iter = LegendrePolynomialIterator(0.5, 4);

julia> collect(iter)
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
  1.0
  0.5
 -0.125
 -0.4375
 -0.2890625

julia> iter = LegendrePolynomialIterator(0.5);

julia> collect(Iterators.take(iter, 5)) # evaluete 5 elements (l = 0:4)
5-element Array{Float64,1}:
  1.0
  0.5
 -0.125
 -0.4375
 -0.2890625

julia> collect(Iterators.take(Iterators.drop(iter, 100), 5)) # evaluate Pl for l = 100:104
5-element Array{Float64,1}:
 -0.0605180259618612
  0.02196749072249231
  0.08178451892628381
  0.05963329258495025
 -0.021651535258316177
```
"""
struct LegendrePolynomialIterator{T, L <: Union{Integer, Nothing}, V}
	x :: V
	lmax :: L
	function LegendrePolynomialIterator{T,L,V}(x::V, lmax::L) where {T, L <: Union{Integer, Nothing}, V}
		checkdomain(x)
		new{T,L,V}(x, lmax)
	end
end
LegendrePolynomialIterator(x) = LegendrePolynomialIterator{polytype(x), Nothing, typeof(x)}(x, nothing)
function LegendrePolynomialIterator(x, lmax)
	assertnonnegative(lmax)
	LegendrePolynomialIterator{polytype(x), typeof(lmax), typeof(x)}(x, lmax)
end

Base.eltype(::Type{<:LegendrePolynomialIterator{T}}) where {T} = T
Base.IteratorSize(::Type{<:LegendrePolynomialIterator{<:Any,Nothing}}) = Base.IsInfinite()

function Base.iterate(::LegendrePolynomialIterator{T}) where {T}
	return one(T), (1, (zero(T), one(T)))
end
function Base.iterate(iter::LegendrePolynomialIterator{T,Nothing}, state) where {T}
	l, (Plm2, Plm1) = state
	x = iter.x
	Pl = Pl_recursion(T, l, Plm1, Plm2, x)
	return Pl, (l+1, (Plm1, Pl))
end
function Base.iterate(iter::LegendrePolynomialIterator{T}, state) where {T}
	l, (Plm2, Plm1) = state
	l > iter.lmax && return nothing
	x = iter.x
	Pl = Pl_recursion(T, l, Plm1, Plm2, x)
	return Pl, (l+1, (Plm1, Pl))
end

Base.length(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = iter.lmax + 1
Base.size(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = (length(iter),)
Base.axes(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = (0:iter.lmax, )
Base.keys(iter::LegendrePolynomialIterator{<:Any,<:Integer}) = 0:iter.lmax

function Base.collect(iter::LegendrePolynomialIterator{T,<:Integer}) where {T}
	arr = zeros(T, 0:iter.lmax)
	collectPl!(arr, iter.x; lmax = iter.lmax)
end

Base.copy(iter::LegendrePolynomialIterator) = typeof(iter)(iter.x, iter.lmax)

"""
	Pl(x, l::Integer)

Compute the Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and the degree `l`

# Examples
```jldoctest
julia> Pl(1, 2)
1.0

julia> Pl(0.5, 20)
-0.04835838106737356
```
"""
function Pl(x, l::Integer)
	iter = LegendrePolynomialIterator(x)
	d = Iterators.drop(iter, l) # 0 to l-1
	first(d)
end

function doublefactorial(T, n)
	p = convert(T, one(n))
	for i in n:-2:1
		p *= convert(T, i)
	end
	convert(T, p)
end

function _checkvalues_m(x, l, m)
	assertnonnegative(l)
	checkdomain(x)
	m >= 0 || throw(ArgumentError("coefficient m must be >= 0"))
end

Base.@propagate_inbounds function _unsafePlm!(cache, x, l, m)
	# unsafe, assumes 1-based indexing
	checklength(cache, l - m + 1)
	# handle m == 0
	collectPl!(cache, x, lmax = l - m)

	# handle m > 0
	for mi in 1:m
		# We denote the terms as P_mi_li
	
		# li == mi
		P_mim1_mim1 = cache[1]
		P_mi_mi = Plm_recursion_m(eltype(cache), mi, mi, P_mim1_mim1, x)
		cache[1] = P_mi_mi

		if l > m
			# li == mi + 1
			# P_mim1_mi = cache[2]
			# P_mi_mip1 = Plm_recursion(eltype(cache), mi + 1, mi, P_mi_mi, P_mim1_mi, nothing, x)
			# cache[2] = P_mi_mip1
			P_mi_lim2 = zero(x)
			P_mi_lim1 = cache[1]
			P_mi_li = Plm_recursion(eltype(cache), mi + 1, mi, P_mi_lim1, P_mi_lim2, x)
			cache[2] = P_mi_li

			for li in mi+2:min(l, l - m + mi)
				P_mi_lim2 = cache[li - mi - 1]
				P_mi_lim1 = cache[li - mi]
				P_mi_li = Plm_recursion(eltype(cache), li, mi, P_mi_lim1, P_mi_lim2, x)
				cache[li - mi + 1] = P_mi_li
			end
		end
	end
	nothing
end

"""
	Plm(x, l::Integer, m::Integer, [cache::AbstractVector])

Compute the associatedLegendre polynomial ``P_\\ell,m(x)``.
Optionally a pre-allocated vector `cache` may be provided, which must have a minimum length of `l - m + 1` 
and may be overwritten during the computation.

The coefficient `m` must be non-negative. For `m == 0` this function just returns 
Legendre polynomials.

# Examples

```jldoctest
julia> Plm(0.5, 3, 2)
5.625

julia> Plm(0.5, 4, 0) == Pl(0.5, 4)
true
```
"""
Base.@propagate_inbounds function Plm(x, l::Integer, m::Integer, 
	A = begin
		_checkvalues_m(x, l, m)
		# do not allocate A if the value is trivially zero
		if l < m
			return zero(polytype(x))
		end
		zeros(polytype(x), l - m + 1)
	end
	)

	_checkvalues_m(x, l, m)
	# check if the value is trivially zero in case A is provided in the function call
	if l < m
		return zero(eltype(A))
	end

	cache = OffsetArrays.no_offset_view(A)
	# function barrier, as no_offset_view may be type-unstable
	_unsafePlm!(cache, x, l, m)

	return cache[l - m + 1]
end

function _checkvalues(x, l, n)
	assertnonnegative(l)
	checkdomain(x)
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
		_checkvalues(x, l, n)
		# do not allocate A if the value is trivially zero
		if l < n
			return zero(polytype(x))
		end
		zeros(polytype(x), l - n + 1)
	end
	)

	_checkvalues(x, l, n)
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
	collectPl!(v::AbstractVector, x; [lmax::Integer = length(v) - 1])

Compute the Legendre Polynomials ``P_\\ell(x)`` for the argument `x` and all degrees `l` in `0:lmax`,
and store them in `v`.

At output `v[firstindex(v) + l] == Pl(x,l)`.

# Examples
```jldoctest
julia> v = zeros(4);

julia> collectPl!(v, 0.5)
4-element Vector{Float64}:
  1.0
  0.5
 -0.125
 -0.4375

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
function collectPl!(v::AbstractVector, x; lmax::Integer  = length(v) - 1)
	checklength(v, lmax + 1)

	iter = LegendrePolynomialIterator(x, lmax)

	for (l, Pl) in pairs(iter)
		v[l + firstindex(v)] = Pl
	end

	v
end

"""
	collectPl(x; lmax::Integer)

Compute the Legendre Polynomial ``P_\\ell(x)`` for the argument `x` and all degrees `l` in `0:lmax`.
Return `P` with indices `0:lmax`, where `P[l] == Pl(x,l)`

# Examples
```jldoctest
julia> collectPl(0.5, lmax = 4)
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
  1.0
  0.5
 -0.125
 -0.4375
 -0.2890625
```
"""
collectPl(x; lmax::Integer) = collect(LegendrePolynomialIterator(x, lmax))

"""
	collectPlm(x; lmax::Integer, m::Integer)

Compute the associated Legendre Polynomial ``P_\\ell,m(x)`` for the argument `x` and all degrees `l = 0:lmax`. 

The coefficient `m` must be greater than or equal to zero.

Returns `v` with indices `0:lmax`, where `v[l] == Plm(x, l, m)`.

# Examples

```jldoctest
julia> collectPlm(0.5, lmax = 3, m = 2)
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
 0.0
 0.0
 2.25
 5.625
```
"""
function collectPlm(x; lmax::Integer, m::Integer)
	assertnonnegative(lmax)
	m >= 0 || throw(ArgumentError("coefficient m must be >= 0"))
	v = zeros(polytype(x), 0:lmax)
	if lmax >= m
		collectPlm!(parent(v), x; lmax = lmax, m = m)
	end
	v
end

"""
	collectPlm!(v::AbstractVector, x; lmax::Integer, m::Integer)

Compute the associated Legendre Polynomial ``P_\\ell,m(x)`` for the argument `x` and all degrees `l = 0:lmax`, 
and store the result in `v`. 

The coefficient `m` must be greater than or equal to zero.

At output, `v[l + firstindex(v)] == Plm(x, l, m)` for `l = 0:lmax`.

# Examples

```jldoctest
julia> v = zeros(4);

julia> collectPlm!(v, 0.5, lmax = 3, m = 2)
4-element Vector{Float64}:
 0.0
 0.0
 2.25
 5.625
```
"""
function collectPlm!(v, x; lmax::Integer, m::Integer)
	assertnonnegative(lmax)
	m >= 0 || throw(ArgumentError("coefficient m must be >= 0"))
	checklength(v, lmax + 1)
	
	# trivially zero for l < m
	fill!((@view v[(0:m-1) .+ firstindex(v)]), zero(eltype(v)))
	# populate the other elements
	@inbounds Plm(x, lmax, m, @view v[(m:lmax) .+ firstindex(v)])
	
	v
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

``jldoctest
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
