module LegendrePolynomials

using HyperDualNumbers
using OffsetArrays

export Pl,
	Pl!,
	dPl,
	dPl!,
	d2Pl,
	d2Pl!,
	Pl_dPl,
	Pl_dPl!,
	Pl_d2Pl,
	Pl_d2Pl!,
	dPl_d2Pl,
	dPl_d2Pl!,
	Pl_dPl_d2Pl,
	Pl_dPl_d2Pl!

##########################################################################################
# Legendre Polynomials and derivatives up to order 2

@inline function Pl_recursion(ℓ, Plm1, Plm2, x)
	((2*(ℓ-1)+1)* x * Plm1 - (ℓ-1) * Plm2)/ℓ
end

function Pl_hyperdual_allmodes(x::Real, lmax::Integer)
	T = typeof(float(x))
	arr = `OffsetArray`{Hyper{T}}(undef,0:lmax)
	Pl_hyperdual_allmodes!(arr,x,lmax)
	return arr
end

function Pl_hyperdual_allmodes!(arr::OffsetVector{Hyper{TR},Vector{Hyper{TR}}}, x::Real,
	lmax::Integer = maximum(axes(arr,1))) where {TR<:Real}

	# Check if x is within limits
	abs(x) > 1 && throw(DomainError(x,
		"Legendre Polynomials are defined for arguments lying in -1 ⩽ x ⩽ 1"))

	if (lmax≥0)
		arr[0] = hyper(one(x))
	end
	if (lmax≥1)
		hdx = hyper(x,one(x),one(x),zero(x))
		arr[1] = hdx
	end
	@inbounds for ℓ in 2:lmax
		arr[ℓ] = Pl_recursion(ℓ,arr[ℓ-1],arr[ℓ-2],hdx)
	end
	arr
end

# Compute Pl without storing intermediate values
function Pl_hyperdual(x::Real, l::Integer)
	# Check if x is within limits
	abs(x) > 1 && throw(DomainError(x,
		"Legendre Polynomials are defined for arguments lying in -1 ⩽ x ⩽ 1"))

	hd1 = hyper(one(x))
	hdx = hyper(x,one(x),one(x),zero(x))

	if l==0
		return hd1
	elseif l==1
		return hdx
	end

	Plm1 = hdx; Plm2 = hd1
	Pl = Plm1

	for ℓ in 2:l
		Pl = Pl_recursion(ℓ,Plm1,Plm2,hdx)
		Plm2,Plm1 = Plm1,Pl
	end

	Pl	
end

function Pl_derivatives_allmodes!(arr::AbstractVector{<:Real}, x::Real,
	lmax::Integer = maximum(axes(arr,1)), deriv=0)

	# Only compute Pl
	arr_hyperdual = view(Pl_hyperdual_allmodes(x,lmax),0:lmax)

	@inbounds @. arr[0:lmax] = realpart(arr_hyperdual)
	arr
end

function Pl_derivatives_allmodes!(arr::AbstractMatrix{<:Real}, x::Real,
	lmax::Integer = maximum(axes(arr,1)),
	deriv::Integer = maximum(axes(arr,2)))

	arr_hyperdual = view(Pl_hyperdual_allmodes(x,lmax),0:lmax)
	
	@inbounds for (d_order,f) in zip(0:deriv,(realpart,eps1,eps1eps2))
		@. arr[0:lmax,d_order] = f(arr_hyperdual)
	end
	arr
end

function Pl_derivatives_allmodes(x::Real, lmax::Integer, deriv::Integer = 0)
	arr = zeros(0:lmax, 0:deriv)
	Pl_derivatives_allmodes!(arr, x, lmax, deriv)
	arr
end

# Convenience functions

"""
	Pl(x::Real; lmax::Integer)

Computes the Legendre Polynomials ``P_l(x)`` for the argument `x` and `l = 0:lmax`.
Returns an `OffsetArray` `P` with indices `0:lmax`, where `P[l] == Pl(x,l)`
"""
Pl(x::Real; lmax::Integer) = Pl_derivatives_allmodes(x,lmax,0)[:,0]

"""
	Pl(x::Real, l::Integer)

Computes the Legendre Polynomials ``P_l(x)`` for the argument `x` and the degree `l`
"""
function Pl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x, l)
	realpart(Ph)
end

"""
	dPl(x::Real; lmax::Integer)

Computes the first derivatives of Legendre Polynomials ``d_xP_l(x)`` for the 
argument `x` and `l = 0:lmax`.
Returns an `OffsetArray` `dP` with indices `0:lmax`, where `dP[l] = dPl(x,l)`
"""
dPl(x::Real; lmax::Integer) = Pl_derivatives_allmodes(x,lmax,1)[:,1]

"""
	dPl(x::Real, l::Integer)

Computes the first derivatives of Legendre Polynomials ``d_xP_l(x)`` for the 
argument `x` and the degree `l`
"""
function dPl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x, l)
	eps1(Ph)
end

"""
	d2Pl(x::Real; lmax::Integer)

Computes the second derivatives of Legendre Polynomials ``d_x^2P_l(x)`` for the 
argument `x` and `l = 0:lmax`
Returns an `OffsetArray` `d2P` with indices `0:lmax`, where `d2P[l] == d2Pl(x,l)`
"""
d2Pl(x::Real; lmax::Integer) = Pl_derivatives_allmodes(x, lmax, 2)[:,2]

"""
	d2Pl(x::Real, l::Integer)

Computes the second derivatives of Legendre Polynomials ``d_x^2P_l(x)`` for the 
argument `x` and the degree `l`
"""
function d2Pl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x, l)
	eps1eps2(Ph)
end

"""
	Pl_dPl(x::Real; lmax::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their derivatives ``d_xP_l(x)``
for the argument `x` and `l = 0:lmax`
Returns `OffsetArray`s `P` and `dP` with indices `0:lmax`, where `P[l] == Pl(x,l)` and 
`dP[l] == dPl(x,l)`.
"""
function Pl_dPl(x::Real; lmax::Integer)
	P = Pl_derivatives_allmodes(x,lmax,1)
	P[:,0], P[:,1]
end

"""
	Pl_dPl(x::Real, l::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their first derivatives 
``d_xP_l(x)`` for the argument `x` and the degree `l`.
"""
function Pl_dPl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph), eps1(Ph)
end

"""
	Pl_dPl_d2Pl(x::Real; lmax::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their first derivatives 
``d_xP_l(x)`` and second derivatives ``d_x^2P_l(x)`` for the argument `x` and `l = 0:lmax`.
Returns `OffsetArray`s `P`, `dP` and `d2P` with indices `0:lmax`, 
where `P[l] == Pl(x,l)`, `dP[l] == dPl(x,l)` and `d2P[l] == d2Pl(x,l)`
"""
function Pl_dPl_d2Pl(x::Real; lmax::Integer)
	P = Pl_derivatives_allmodes(x, lmax, 2)
	P[:,0], P[:,1], P[:,2]
end

"""
	Pl_dPl_d2Pl(x::Real, l::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their first derivatives 
``d_xP_l(x)`` and second derivatives ``d_x^2P_l(x)`` for the argument ``x`` and the degree ``l``
"""
function Pl_dPl_d2Pl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph), eps1(Ph), eps1eps2(Ph)
end

"""
	dPl_d2Pl(x::Real; lmax::Integer)

Computes the the first and second derivatives of Legendre Polynomials ``d_xP_l(x)``
and ``d_x^2P_l(x)`` for the argument `x` and `l = 0:lmax`

Returns `OffsetArray`s `dP` and `d2P` with indices `0:lmax`, where `dP[l] == dPl(x,l)` and 
`d2P[l] == d2Pl(x,l)`
"""
function dPl_d2Pl(x::Real; lmax::Integer)
	P = Pl_derivatives_allmodes(x, lmax, 2)
	P[:,1], P[:,2]
end

"""
	dPl_d2Pl(x::Real, l::Integer)

Computes the the first and second derivatives of Legendre Polynomials ``d_xP_l(x)``
and ``d_x^2P_l(x)`` for the argument `x` and the degree `l`
"""
function dPl_d2Pl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x,l)
	eps1(Ph), eps1eps2(Ph)
end

"""
	Pl_d2Pl(x::Real; lmax::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their second derivatives 
``d_x^2P_l(x)`` for the argument `x` and `l = 0:lmax`

Returns `OffsetArray`s `P` and `d2P` with indices `0:lmax`, where `P[l] == Pl(x,l)` and 
`d2P[l] == d2Pl(x,l)`
"""
function Pl_d2Pl(x::Real; lmax::Integer)
	P = Pl_derivatives_allmodes(x, lmax, 2)
	P[:,0], P[:,2]
end

"""
	Pl_d2Pl(x::Real, l::Integer)

Computes the the Legendre Polynomials ``P_l(x)`` and their second derivatives 
``d_x^2P_l(x)`` for the argument `x` and the degree `l`
"""
function Pl_d2Pl(x::Real, l::Integer)
	Ph = Pl_hyperdual(x, l)
	realpart(Ph), eps1eps2(Ph)
end

"""
	Pl!(arr::AbstractArray, x::Real; [lmax = maximum(axes(arr,1))])

Computes the Legendre Polynomials ``P_l(x)`` for the argument `x` and `l = 0:lmax`,
and saves it in `arr`. Assumes that `arr` has `0`-based indexing.

All dimensions of `arr` should have indices starting from `0`.
The axes of `arr` should be `(0:l, 0:n)` where `l >= lmax` and `n >= 0`

At output, `arr[l,0] == Pl(x,l)`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
Pl!(arr::AbstractVecOrMat{<:Real}, x::Real; lmax::Integer = maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr,x,lmax,0)

"""
	Pl_dPl!(arr, x; [lmax = maximum(axes(arr,1))])

Computes the Legendre Polynomials ``P_l(x)`` and their derivatives ``d_xP_l(x)``
for the argument `x` and `l = 0:lmax`, and saves them in `arr`. 
Assumes that `arr` has `0`-based indexing.

All dimensions of `arr` should have indices starting from `0` 
The axes of `arr` should be `(0:l, 0:n)` where `l >= lmax` and ``n >= 1``

At output, `arr[l,0] == Pl(x,l)` and `arr[l,1] == dPl(x,l)`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
Pl_dPl!(arr::AbstractMatrix{<:Real}, x::Real; lmax::Integer = maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr, x, lmax, min(1,maximum(axes(arr,2))))

"""
	dPl!(arr::AbstractArray, x::Real; [lmax = maximum(axes(arr,1))])

Computes the first derivatives of Legendre Polynomials ``d_xP_l(x)``
for the argument `x` and `l = 0:lmax`, and saves them in `arr`. 
Assumes that arr has `0`-based indexing.

The first dimension of arr should be `0:l`, where `l >= lmax`
The second dimension -- if `arr` is 2D -- should contain the index `1`

At output, `arr[l] == dPl(x,l)` if `arr` is a `Vector`
or `arr[l,1] == dPl(x,l)` if `arr` is a `Matrix`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
function dPl!(arr::AbstractVector{<:Real}, x::Real; lmax::Integer = maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:1)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)
	@views @. arr = P[:,1]
	arr
end

function dPl!(arr::AbstractMatrix{<:Real}, x::Real;
	lmax::Integer = maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:1)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)
	@views @. arr[:,1] = P[:,1]
	arr
end

"""
	Pl_dPl_d2Pl!(arr::AbstractMatrix, x::Real; [lmax = maximum(axes(arr,1))])

Computes the Legendre Polynomials ``P_l(x)`` and their first and second 
derivatives ``d_xP_l(x)`` and ``d_x^2P_l(x)``, for the argument `x` and `l = 0:lmax`, 
and saves them in `arr`. Assumes that `arr` has `0`-based indexing.

All dimensions of `arr` should have indices starting from `0`
The axes of `arr` should be `(0:l, 0:n)` where `l >= lmax` and `n >= 2`

At output, `arr[l,0] == Pl(x,l)`, `arr[l,1] == dPl(x,l)` and `arr[l,2] == d2Pl(x,l)`

The optional keyword argument `lmax` can specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
Pl_dPl_d2Pl!(arr::AbstractMatrix{<:Real}, x::Real; lmax::Integer = maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr, x, lmax, min(2,maximum(axes(arr,2))))

"""
	d2Pl!(arr::AbstractArray, x::Real; [lmax = maximum(axes(arr,1))])

Computes the second derivatives of Legendre Polynomials ``d_x^2P_l(x)``
for the argument `x` and `l = 0:lmax`, and saves them in `arr`. 
Assumes that `arr` has `0`-based indexing.

The first dimension of `arr` should be `0:l`, where `l >= lmax`
The second dimension -- if `arr` is 2D -- should contain the index `2`

At output, `arr[l] == d2Pl(x,l)` if `arr` is a `Vector`
or `arr[l,2] == d2Pl(x,l)` if `arr` is a `Matrix`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
function d2Pl!(arr::AbstractVector{<:Real}, x::Real; lmax::Integer = maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)
	@views @. arr = P[:,2]
	arr
end

function d2Pl!(arr::AbstractMatrix{<:Real},
	x::Real; lmax::Integer = maximum(axes(arr,1)))

	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)
	@views @. arr[:,2] = P[:,2]
	arr
end

"""
	dPl_d2Pl!(arr::AbstractMatrix, x::Real; [lmax = maximum(axes(arr,1))])

Computes the first and second derivatives ``d_xP_l(x)`` and ``d_x^2P_l(x)``, 
for the argument `x` and `l = 0:lmax`, 
and saves them in `arr`. Assumes that `arr` has `0`-based indexing.

The first dimension of arr should be `0:l`, where `l >= lmax`
The second dimension should be `1:n` where `n >= 2`

At output, `arr[l,1] == dPl(x,l)` and `arr[l,2] == d2Pl(x,l)`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
function dPl_d2Pl!(arr::AbstractMatrix{<:Real}, x::Real;
	lmax::Integer = maximum(axes(arr,1)))

	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)

	@views @. arr[:,1:2] = P[:,1:2]
	
	arr
end

"""
	Pl_d2Pl!(arr::AbstractMatrix, x::Real; [lmax = maximum(axes(arr,1))])

Computes the first and second derivatives ``d_xP_l(x)`` and ``d_x^2P_l(x)``, 
for the argument `x` and `l = 0:lmax`, 
and saves them in `arr`. Assumes that `arr` has `0`-based indexing.

The first dimension of arr should be `0:l`, where `l >= lmax`
The second dimension should be `0:n` where `n >= 2`

At output, `arr[l,0] == dPl(x,l)` and `arr[l,2] == d2Pl(x,l)`

The optional keyword argument `lmax` may specify the range of ``l``'s to compute.
It defaults to `lmax = maximum(axes(arr,1))`
"""
function Pl_d2Pl!(arr::AbstractMatrix{<:Real}, x::Real;
	lmax::Integer = maximum(axes(arr,1)))

	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P, x; lmax = lmax)

	@views @. arr[:,0] = P[:,0]
	@views @. arr[:,2] = P[:,2]
	
	arr
end

end