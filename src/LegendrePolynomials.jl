module LegendrePolynomials
using HyperDualNumbers,OffsetArrays

export Pl,Pl!,
dPl,dPl!,
d2Pl,d2Pl!,
Pl_dPl,Pl_dPl!,
Pl_d2Pl,Pl_d2Pl!,
dPl_d2Pl,dPl_d2Pl!,
Pl_dPl_d2Pl, Pl_dPl_d2Pl!

##########################################################################################
# Legendre Polynomials and derivatives up to order 2

@inline function Pl_recursion(ℓ,Plm1,Plm2,x)
	((2*(ℓ-1)+1)* x * Plm1 - (ℓ-1) * Plm2)/ℓ
end

function Pl_hyperdual_allmodes(x::Real,lmax::Integer)
	T = typeof(float(x))
	arr = OffsetArray{Hyper{T}}(undef,0:lmax)
	Pl_hyperdual_allmodes!(arr,x,lmax)
	return arr
end

function Pl_hyperdual_allmodes!(arr::OffsetVector{Hyper{TR},Vector{Hyper{TR}}},x::Real,
	lmax::Integer=maximum(axes(arr,1))) where {TR<:Real}

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
function Pl_hyperdual(x::Real,l::Integer)
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

function Pl_derivatives_allmodes!(arr::AbstractVector{<:Real},x::Real,
	lmax::Integer=maximum(axes(arr,1)),deriv=0)

	# Only compute Pl
	arr_hyperdual = view(Pl_hyperdual_allmodes(x,lmax),0:lmax)

	@inbounds @. arr[0:lmax] = realpart(arr_hyperdual)
	arr
end

function Pl_derivatives_allmodes!(arr::AbstractMatrix{<:Real},x::Real,
	lmax::Integer=maximum(axes(arr,1)),
	deriv::Integer=maximum(axes(arr,2)))

	arr_hyperdual = view(Pl_hyperdual_allmodes(x,lmax),0:lmax)
	
	@inbounds for (d_order,f) in zip(0:deriv,(realpart,eps1,eps1eps2))
		@. arr[0:lmax,d_order] = f(arr_hyperdual)
	end
	arr
end

function Pl_derivatives_allmodes(x::Real,lmax::Integer,deriv::Integer=0)
	arr = zeros(0:lmax,0:deriv)
	Pl_derivatives_allmodes!(arr,x,lmax,deriv)
	arr
end

# Convenience functions

"""
	Pl(x::Real;lmax::Integer)

Computes the Legendre Polynomials Pl(x) for the argument x and l=0:lmax.
Returns an OffsetArray P[0:lmax], where P[l] = Pl(x)
"""
Pl(x::Real;lmax::Integer) = Pl_derivatives_allmodes(x,lmax,0)[:,0]

"""
	Pl(x::Real,l::Integer)

Computes the Legendre Polynomials Pl(x) for the argument x and the degree l
"""
function Pl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph)
end

"""
	dPl(x::Real;lmax::Integer)

Computes the first derivatives of Legendre Polynomials dₓPl(x) for the 
argument x and l=0:lmax.
Returns an OffsetArray dP[0:lmax], where `dP[l]` = dₓPl(x)
"""
dPl(x::Real;lmax::Integer) = Pl_derivatives_allmodes(x,lmax,1)[:,1]

"""
	dPl(x::Real,l::Integer)

Computes the first derivatives of Legendre Polynomials dₓPl(x) for the 
argument x and the degree l
"""
function dPl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	eps1(Ph)
end

"""
	d2Pl(x::Real;lmax::Integer)

Computes the second derivatives of Legendre Polynomials dₓ²Pl(x) for the 
argument x and l=0:lmax
Returns an OffsetArray d2P[0:lmax], where d2P[l] = dₓ²Pl(x)
"""
d2Pl(x::Real;lmax::Integer) = Pl_derivatives_allmodes(x,lmax,2)[:,2]

"""
	d2Pl(x::Real,l::Integer)

Computes the second derivatives of Legendre Polynomials dₓ²Pl(x) for the 
argument x and the degree l
"""
function d2Pl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	eps1eps2(Ph)
end

"""
	Pl_dPl(x::Real;lmax::Integer)

Computes the the Legendre Polynomials Pl(x) and their derivatives dₓPl(x) 
for the argument x and l=0:lmax
Returns OffsetArrays P[0:lmax] and dP[0:lmax], where P[l] = Pl(x) and 
dP[l] = dₓPl(x)
"""
function Pl_dPl(x::Real;lmax::Integer)
	P = Pl_derivatives_allmodes(x,lmax,1)
	P[:,0],P[:,1]
end

"""
	Pl_dPl(x::Real,l::Integer)

Computes the the Legendre Polynomials Pl(x) and their first derivatives 
dPl(x) for the argument x and the degree l
"""
function Pl_dPl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph),eps1(Ph)
end

"""
	Pl_dPl_d2Pl(x::Real;lmax::Integer)

Computes the the Legendre Polynomials Pl(x) and their first derivatives 
dPl(x) and second derivatives dₓ²Pl(x) for the argument x and l=0:lmax
Returns OffsetArrays P[0:lmax], dP[0:lmax] and d2P[0:lmax], 
where P[l] = Pl(x), dP[l] = dₓPl(x) and d2P[l] = dₓ²Pl(x)
"""
function Pl_dPl_d2Pl(x::Real;lmax::Integer)
	P = Pl_derivatives_allmodes(x,lmax,2)
	P[:,0],P[:,1],P[:,2]
end

"""
	Pl_dPl_d2Pl(x::Real,l::Integer)

Computes the the Legendre Polynomials Pl(x) and their first derivatives 
dPl(x) and second derivatives dₓ²Pl(x) for the argument x and the degree l
"""
function Pl_dPl_d2Pl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph),eps1(Ph),eps1eps2(Ph)
end

"""
	dPl_d2Pl(x::Real;lmax::Integer)

Computes the the first and second derivatives of Legendre Polynomials dPl(x) 
and dₓ²Pl(x) for the argument x and l=0:lmax

Returns OffsetArrays dP[0:lmax] and d2P[0:lmax], where dP[l] = dₓPl(x) and 
d2P[l] = dₓ²Pl(x)
"""
function dPl_d2Pl(x::Real;lmax::Integer)
	P = Pl_derivatives_allmodes(x,lmax,2)
	P[:,1],P[:,2]
end

"""
	dPl_d2Pl(x::Real,l::Integer)

Computes the the first and second derivatives of Legendre Polynomials dPl(x) 
and dₓ²Pl(x) for the argument x and the degree l
"""
function dPl_d2Pl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	eps1(Ph),eps1eps2(Ph)
end

"""
	Pl_d2Pl(x::Real;lmax::Integer)

Computes the the Legendre Polynomials Pl(x) and their second derivatives 
dₓ²Pl(x) for the argument x and l=0:lmax

Returns OffsetArrays P[0:lmax] and d2P[0:lmax], where P[l] = Pl(x) and 
d2P[l] = dₓ²Pl(x)
"""
function Pl_d2Pl(x::Real;lmax::Integer)
	P = Pl_derivatives_allmodes(x,lmax,2)
	P[:,0],P[:,2]
end

"""
	Pl_d2Pl(x::Real,l::Integer)

Computes the the Legendre Polynomials Pl(x) and their second derivatives 
dₓ²Pl(x) for the argument x and the degree l
"""
function Pl_d2Pl(x::Real,l::Integer)
	Ph = Pl_hyperdual(x,l)
	realpart(Ph),eps1eps2(Ph)
end

"""
	Pl!(arr,x;[lmax])

Computes the Legendre Polynomials Pl(x) for the argument x and l=0:lmax,
and saves it in the AbstractArray arr. Assumes that arr has 0-based indexing.

All dimensions of arr should have indices starting from 0
The axes of arr should be (0:l,0:n) where l ⫺ lmax and n ⫺ 0

At output, arr[l,0] = Pl(x)

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
Pl!(arr::AbstractVecOrMat{<:Real},x::Real;lmax::Integer=maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr,x,lmax,0)

"""
	Pl_dPl!(arr,x;[lmax])

Computes the Legendre Polynomials Pl(x) and their derivatives dₓPl(x)
for the argument x and l=0:lmax, and saves them in the AbstractArray arr. 
Assumes that arr has 0-based indexing.

All dimensions of arr should have indices starting from 0 
The axes of arr should be (0:l,0:n) where l ⫺ lmax and n ⫺ 1

At output, arr[l,0] = Pl(x) and arr[l,1] = dₓPl(x)

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
Pl_dPl!(arr::AbstractMatrix{<:Real},x::Real;lmax::Integer=maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr,x,lmax,min(1,maximum(axes(arr,2))))

"""
	dPl!(arr,x;[lmax])

Computes the first derivatives of Legendre Polynomials dₓPl(x)
for the argument x and l=0:lmax, and saves them in the AbstractArray arr. 
Assumes that arr has 0-based indexing.

The first dimension of arr should be 0:l, where l ⫺ lmax
The second dimension -- if arr is 2D -- should contain the index 1

At output, arr[l] = dₓPl(x) if arr is a Vector
or arr[l,1] = dₓPl(x) if arr is a Matrix

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
function dPl!(arr::AbstractVector{<:Real},x::Real;lmax::Integer=maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:1)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)
	@inbounds @. arr[:] = P[:,1]
	arr
end

function dPl!(arr::AbstractMatrix{<:Real},x::Real;
	lmax::Integer=maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:1)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)
	@inbounds @. arr[:,1] = P[:,1]
	arr
end

"""
	Pl_dPl_d2Pl!(arr,x;[lmax])

Computes the Legendre Polynomials Pl(x) and their first and second 
derivatives dₓPl(x) and dₓ²Pl(x), for the argument x and l=0:lmax, 
and saves them in the AbstractArray arr. Assumes that arr has 0-based indexing.

All dimensions of arr should have indices starting from 0
The axes of arr should be (0:l,0:n) where l ⫺ lmax and n ⫺ 2

At output, arr[l,0] = Pl(x), arr[l,1] = dₓPl(x) and arr[l,2] = dₓ²Pl(x)

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
Pl_dPl_d2Pl!(arr::AbstractMatrix{<:Real},x::Real;lmax::Integer=maximum(axes(arr,1))) = 
	Pl_derivatives_allmodes!(arr,x,lmax,min(2,maximum(axes(arr,2))))

"""
	d2Pl!(arr,x;[lmax])

Computes the second derivatives of Legendre Polynomials dₓ²Pl(x)
for the argument x and l=0:lmax, and saves them in the AbstractArray arr. 
Assumes that arr has 0-based indexing.

The first dimension of arr should be 0:l, where l ⫺ lmax
The second dimension -- if arr is 2D -- should contain the index 2

At output, arr[l] = dₓ²Pl(x) if arr is a Vector
or arr[l,2] = dₓ²Pl(x) if arr is a Matrix

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
function d2Pl!(arr::AbstractVector{<:Real},x::Real;lmax::Integer=maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)
	@. arr[:] = P[:,2]
	arr
end

function d2Pl!(arr::AbstractMatrix{<:Real},
	x::Real;lmax::Integer=maximum(axes(arr,1)))
	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)
	@. arr[:,2] = P[:,2]
	arr
end

"""
	dPl_d2Pl!(arr,x;[lmax])

Computes the first and second derivatives dₓPl(x) and dₓ²Pl(x), 
for the argument x and l=0:lmax, 
and saves them in the AbstractArray arr. Assumes that arr has 0-based indexing.

The first dimension of arr should be 0:l, where l ⫺ lmax
The second dimension should be (1:n) where n ⫺ 2

At output, arr[l,1] = dₓPl(x) and arr[l,2] = dₓ²Pl(x)

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
function dPl_d2Pl!(arr::AbstractMatrix{<:Real},x::Real;
	lmax::Integer=maximum(axes(arr,1)))

	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)

	@. arr[:,1] = P[:,1]
	@. arr[:,2] = P[:,2]
	
	arr
end

"""
	Pl_d2Pl!(arr,x;[lmax])

Computes the first and second derivatives dₓPl(x) and dₓ²Pl(x), 
for the argument x and l=0:lmax, 
and saves them in the AbstractArray arr. Assumes that arr has 0-based indexing.

The first dimension of arr should be 0:l, where l ⫺ lmax
The second dimension should be (0:n) where n ⫺ 2

At output, arr[l,0] = dₓPl(x) and arr[l,2] = dₓ²Pl(x)

The optional keyword argument lmax can specify the range of l's to compute.
It defaults to lmax = maximum(axes(arr,1))
"""
function Pl_d2Pl!(arr::AbstractMatrix{<:Real},x::Real;
	lmax::Integer=maximum(axes(arr,1)))

	P = zeros(axes(arr,1),0:2)
	Pl_dPl_d2Pl!(P,x;lmax=lmax)

	@. arr[:,0] = P[:,0]
	@. arr[:,2] = P[:,2]
	
	arr
end

end