module LegendrePolynomials

using OffsetArrays

export Pl
export collectPl
export collectPl!

function checkdomain(x)
	abs(x) > 1 && throw(DomainError(x,
		"Legendre Polynomials are defined for arguments lying in -1 ⩽ x ⩽ 1"))
end

function checksize(arr, lmax)
	maximum(axes(arr,1)) >= lmax || throw(ArgumentError("array is not large enough to store all values"))
end

@inline function Pl_recursion(::Type{T}, ℓ, Plm1, Plm2, x) where {T}
	# relation is valid from ℓ = 1
	Pl = ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ
	convert(T, Pl)
end

polytype(x) = typeof(x*x/1)

"""
	LegendrePolynomialIterator(x, [lmax::Integer])

Return an iterator that generates the values of the Legendre polynomials ``P_l(x)`` for the given `x`.
If `lmax` is specified then only the values of ``P_l(x)`` from `0` to `lmax` are returned.

# Examples
```jldoctest
julia> iter = LegendrePolynomialIterator(0.5, 4);

julia> collect(iter)
5-element OffsetArray(::Array{Float64,1}, 0:4) with eltype Float64 with indices 0:4:
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
end
LegendrePolynomialIterator(x) = LegendrePolynomialIterator{polytype(x), Nothing, typeof(x)}(x, nothing)
function LegendrePolynomialIterator(x, lmax)
	lmax >= 0 || throw(ArgumentError("degree must be >= 0"))
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

Compute the Legendre Polynomial ``P_l(x)`` for the argument `x` and the degree `l`

# Examples
```jldoctest
julia> Pl(1, 2)
1.0

julia> Pl(0.5, 20)
-0.04835838106737356
```
"""
function Pl(x, l::Integer)
	# Check if x is within limits
	checkdomain(x)

	iter = LegendrePolynomialIterator(x)
	d = Iterators.drop(iter, l) # 0 to l-1
	first(d)
end

"""
	collectPl!(v::AbstractVector, x; [lmax::Integer = length(v) - 1])

Compute the Legendre Polynomials ``P_l(x)`` for the argument `x` and all degrees `l` in `0:lmax`, 
and store them in `v`.

At output `v[firstindex(v) + l] == Pl(x,l)`.

# Examples
```jldoctest
julia> v = zeros(4);

julia> collectPl!(v, 0.5)
4-element Array{Float64,1}:
  1.0
  0.5
 -0.125
 -0.4375

julia> v = zeros(0:4);

julia> collectPl!(v, 0.5, lmax = 3) # only l from 0 to 3 are populated
5-element OffsetArray(::Array{Float64,1}, 0:4) with eltype Float64 with indices 0:4:
  1.0
  0.5
 -0.125
 -0.4375
  0.0
```
"""
function collectPl!(v::AbstractVector, x; lmax::Integer  = length(v) - 1)
	v_0based = OffsetArray(v, OffsetArrays.Origin(0))
	checksize(v_0based, lmax)

	iter = LegendrePolynomialIterator(x, lmax)
	
	for (l, Pl) in pairs(iter)
		v_0based[l] = Pl
	end

	v
end

"""
	collectPl(x; lmax::Integer)

Compute the Legendre Polynomial ``P_l(x)`` for the argument `x` and all degrees `l` in `0:lmax`.
Return an `OffsetArray` `P` with indices `0:lmax`, where `P[l] == Pl(x,l)`

# Examples
```jldoctest
julia> collectPl(0.5, lmax = 4)
5-element OffsetArray(::Array{Float64,1}, 0:4) with eltype Float64 with indices 0:4:
  1.0
  0.5
 -0.125
 -0.4375
 -0.2890625
```
"""
collectPl(x; lmax::Integer) = collect(LegendrePolynomialIterator(x, lmax))

end