module Legendre
using HyperDualNumbers,OffsetArrays

export Pl,Pl!,Pl_dPl,Pl_dPl!, Pl_dPl_d2Pl, Pl_dPl_d2Pl!
export dPl,d2Pl
export Plm,Plm_cosθ,Plm!,Plm_cosθ!
export Ylm,Ylm!
##########################################################################################
# Legendre Polynomians and derivatives up to order 2

function Pl_hyperdual(x,ℓmax)
	arr = OffsetArray{Hyper{Float64}}(undef,0:ℓmax)
	Pl_hyperdual!(arr,x,ℓmax)
	return arr
end

function Pl_hyperdual!(arr,x,ℓmax=axes(arr,1)[end])
	hd2 = hyper(x,1.0,1.0,0.0)
	if (ℓmax≥0) 
		arr[0] = hyper(1)
	end
	if (ℓmax≥1) 
		arr[1] = hd2
	end
	for ℓ in 2:ℓmax
		arr[ℓ] = (2*(ℓ-1)+1)/ℓ * hd2 * arr[ℓ-1] - (ℓ-1)/ℓ * arr[ℓ-2]
	end
end

function Pl_derivatives!(arr,x,ℓmax=axes(arr,1)[end],deriv=axes(arr,2)[end])
	arr_hyperdual = Pl_hyperdual(x,ℓmax)
	
	for (d_order,f) in zip(0:deriv,(realpart,eps1,eps1eps2))
		arr[0:ℓmax,d_order] .= f.(parent(arr_hyperdual[0:ℓmax]))
	end

end

function Pl_derivatives(x,ℓmax,deriv=0)
	arr = OffsetArray{Float64}(undef,0:ℓmax,0:deriv)
	Pl_derivatives!(arr,x,ℓmax,deriv)
	return arr
end

# Convenience functions

Pl(x;ℓmax) = Pl_derivatives(x,ℓmax,0)[:,0]
Pl_dPl(x;ℓmax) = Pl_derivatives(x,ℓmax,1)
Pl_dPl_d2Pl(x;ℓmax) = Pl_derivatives(x,ℓmax,2)

Pl!(arr,x;ℓmax=axes(arr,1)[end]) = Pl_derivatives!(arr,x,ℓmax,0)
Pl_dPl!(arr,x;ℓmax=axes(arr,1)[end]) = Pl_derivatives!(arr,x,ℓmax,1)
Pl_dPl_d2Pl!(arr,x;ℓmax=axes(arr,1)[end]) = Pl_derivatives!(arr,x,ℓmax,2)

dPl(x;ℓmax) = Pl_derivatives(x,ℓmax,1)[:,1]
d2Pl(x;ℓmax) = Pl_derivatives(x,ℓmax,2)[:,2]

############################################################################################################
############################################################################################################

# Associated Legendre for -2 ≤ m ≤ 2

############################################################################################################
############################################################################################################

function Plm!(arr,x;ℓmax=axes(arr,1)[end])
	dn_Pl = OffsetArray{Float64}(undef,0:ℓmax,0:2)
	Pl_derivatives!(dn_Pl,x,ℓmax,2)

	prefactor = zeros(2)

	for ℓ in 0:ℓmax
		
		prefactor[1] = 1/((ℓ+1)*ℓ)
		prefactor[2] = 1/((ℓ+2)*(ℓ+1)*ℓ*(ℓ-1))

		arr[ℓ,0] = dn_Pl[ℓ,0]

		for m in 1:min(2,ℓ)
			
			arr[ℓ,m] = (-1)^m * (1-x^2)^(m/2) * dn_Pl[ℓ,m]

			arr[ℓ,-m] = (-1)^m * prefactor[m] * arr[ℓ,m]
		end

		for m in -2:-(ℓ+1)
			arr[ℓ,m] = 0.
		end

		for m in (ℓ+1):2
			arr[ℓ,m] = 0.
		end
		
	end
end

function Plm_cosθ!(arr,θ;ℓmax=axes(arr,1)[end])
	dn_Pl = OffsetArray{Float64}(undef,0:ℓmax,0:2)
	Pl_derivatives!(dn_Pl,cos(θ),ℓmax,2)

	prefactor = zeros(2)

	for ℓ in 0:ℓmax
		
		prefactor[1] = 1/((ℓ+1)*ℓ)
		prefactor[2] = 1/((ℓ+2)*(ℓ+1)*ℓ*(ℓ-1))

		arr[ℓ,0] = dn_Pl[ℓ,0]

		for m in 1:min(2,ℓ)
			
			arr[ℓ,m] = (-1)^m * abs(sin(θ))^m * dn_Pl[ℓ,m]

			arr[ℓ,-m] = (-1)^m * prefactor[m] * arr[ℓ,m]
		end

		for m in -2:-(ℓ+1)
			arr[ℓ,m] = 0.
		end

		for m in (ℓ+1):2
			arr[ℓ,m] = 0.
		end
		
	end
end

function Plm(x;ℓmax)
	arr = OffsetArray{Float64}(undef,0:ℓmax,-2:2)
	Plm!(arr,x,ℓmax=ℓmax)
	return arr
end

function Plm_cosθ(θ;ℓmax)
	arr = OffsetArray{Float64}(undef,0:ℓmax,-2:2)
	Plm_cosθ!(arr,θ,ℓmax=ℓmax)
	return arr
end

############################################################################################################
############################################################################################################

# Spherical Harmonics for -2 ≤ m ≤ 2

############################################################################################################
############################################################################################################

function SHnorm(ℓ,m)
	m==0 && return √((2ℓ+1)/4π)
	m==1 && return √((2ℓ+1)/4π/(ℓ*(ℓ+1)))
	m==2 && return √((2ℓ+1)/4π/((ℓ-1)*ℓ*(ℓ+1)*(ℓ+2)))
	m==-1 && return √((2ℓ+1)/4π*(ℓ*(ℓ+1)))
	m==2 && return √((2ℓ+1)/4π*((ℓ-1)*ℓ*(ℓ+1)*(ℓ+2)))
end

function Ylm!(arr,θ,ϕ,ℓmax=axes(arr,1)[end])
	Plm_cosθ!(arr,θ,ℓmax)
	for m in -2:2, ℓ in abs(m):ℓmax
		arr[ℓ,m] *= cis(m*ϕ)*SHnorm(ℓ,m)
	end
end

function Ylm(θ,ϕ,ℓmax)
	arr = OffsetArray{ComplexF64}(undef,0:ℓmax,-2:2)
	Ylm!(arr,θ,ϕ,ℓmax)
	return arr
end


end