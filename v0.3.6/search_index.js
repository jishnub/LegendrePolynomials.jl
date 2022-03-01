var documenterSearchIndex = {"docs":
[{"location":"derivatives/#Derivatives-of-Legendre-Polynomials","page":"Derivatives","title":"Derivatives of Legendre Polynomials","text":"","category":"section"},{"location":"derivatives/#Analytical-recursive-approach","page":"Derivatives","title":"Analytical recursive approach","text":"","category":"section"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Bonnet's recursion formula","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"P_ell(x) = left((2ell-1) x P_ell-1(x) - (ell-1)P_ell - 2(x)right)ell","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"may be differentiated an arbitrary number of times analytically to obtain recursion relations for higher derivatives:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"fracd^n P_ell(x)dx^n = frac(2ell-1)ell left(x fracd^n P_ell-1(x)dx^n +\nn fracd^(n-1) P_ell-1(x)dx^(n-1) right) - frac(ell-1)ell fracd^n P_ell-2(x)dx^n","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"This provides a simultaneous recursion relation in ell as well as n, solving which we may obtain derivatives up to any order. This is the approach used in this package to compute the derivatives of Legendre polynomials.","category":"page"},{"location":"derivatives/#Automatic-diferentiation","page":"Derivatives","title":"Automatic diferentiation","text":"","category":"section"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Julia automatic differentiation framework may be used to compute the derivatives of Legendre polynomials alongside their values. Since the defintions of the polynomials are completely general, they may be called with dual or hyperdual numbers as arguments to evaluate derivarives in one go. We demonstrate one example of this using the package HyperDualNumbers.jl v4:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"DocTestSetup = quote\n\tusing LegendrePolynomials\n  using HyperDualNumbers\nend","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> x = 0.5;\n\njulia> Pl(x, 3)\n-0.4375\n\njulia> using HyperDualNumbers\n\njulia> xh = Hyper(x, one(x), one(x), zero(x));\n\njulia> p = Pl(xh, 3)\n-0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Legendre polynomial P_ell(x) may be obtained using","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> realpart(p)\n-0.4375","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The first derivative dP_ell(x)dx may be obtained as","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁part(p)\n0.375","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The second derivative d^2P_ell(x)dx^2 may be obtained using","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁ε₂part(p)\n7.5","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Something similar may also be evaluated for all l iteratively. For example, the function collectPl evaluates Legendre polynomials for the HyperDualNumber argument for a range of l:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> collectPl(xh, lmax=4)\n5-element OffsetArray(::Vector{Hyper{Float64}}, 0:4) with eltype Hyper{Float64} with indices 0:4:\n                1.0 + 0.0ε₁ + 0.0ε₂ + 0.0ε₁ε₂\n                0.5 + 1.0ε₁ + 1.0ε₂ + 0.0ε₁ε₂\n             -0.125 + 1.5ε₁ + 1.5ε₂ + 3.0ε₁ε₂\n        -0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂\n -0.2890625 - 1.5625ε₁ - 1.5625ε₂ + 5.625ε₁ε₂","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"We may extract the first derivatives by broadcasting the function ε₁part on the array as:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁part.(collectPl(xh, lmax=4))\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  0.0\n  1.0\n  1.5\n  0.375\n -1.5625","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Similarly the function ε₁ε₂part may be used to obtain the second derivatives.","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Several convenience functions to compute the derivatives of Legendre polynomials were available in LegendrePolynomials v0.2, but have been removed in v0.3. The users are encouraged to implement convenience functions to extract the derivatives as necessary. As an exmaple, we may compute the polynomials and their first and second derivatives together as","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> using HyperDualNumbers\n\njulia> function Pl_dPl_d2Pl(x; lmax)\n           xh = Hyper(x, one(x), one(x), zero(x))\n           p = collectPl(xh, lmax = lmax)\n           realpart.(p), ε₁part.(p), ε₁ε₂part.(p)\n       end\nPl_dPl_d2Pl (generic function with 1 method)\n\njulia> Pl_dPl_d2Pl(0.5, lmax = 3)\n([1.0, 0.5, -0.125, -0.4375], [0.0, 1.0, 1.5, 0.375], [0.0, 0.0, 3.0, 7.5])","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = LegendrePolynomials\nDocTestSetup = quote\n\tusing LegendrePolynomials\nend","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Compute Legendre polynomials, Associated Legendre polynomials, and their derivatives using a 3-term recursion relation (Bonnet’s recursion formula).","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"P_ell(x) = left((2ell-1) x P_ell-1(x) - (ell-1)P_ell - 2(x)right)ell","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Currently this package evaluates the standard polynomials that satisfy P_ell(1) = 1 and P_0(x) = 1. These are normalized as","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"int_-1^1 P_m(x) P_n(x) dx = frac22n+1 delta_mn","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"There are six main functions: ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pl(x,l): this evaluates the Legendre polynomial for a given degree l at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectPl(x; lmax): this evaluates all the polynomials for l lying in 0:lmax at the argument x. As before the argument needs to lie in the domain of validity. Functionally this is equivalent to Pl.(x, 0:lmax), except collectPl evaluates the result in one pass, and is therefore faster. There is also the in-place version collectPl! that uses a pre-allocated array.\nPlm(x, l, m): this evaluates the associated Legendre polynomial P_ellm(x) at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectPlm(x; m, lmax): this evaluates the associated Legendre polynomials with coefficient m for l = 0:lmax. There is also an in-place version collectPlm! that uses a pre-allocated array.\ndnPl(x, l, n): this evaluates the n-th derivative of the Legendre polynomial P_ell(x) at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectdnPl(x; n, lmax): this evaluates the n-th derivative of all the Legendre polynomials for l = 0:lmax. There is also an in-place version collectdnPl! that uses a pre-allocated array.","category":"page"},{"location":"#Quick-Start","page":"Introduction","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the Legendre polynomial for one l at an argumentx as Pl(x, l):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> Pl(0.5, 3)\n-0.4375","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the associated Legendre Polynomial one l,m pair as Plm(x, l, m):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> Plm(0.5, 3, 2)\n5.625","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the nth derivative for one l as dnPl(x, l, n):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> dnPl(0.5, 3, 2)\n7.5","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate all the polynomials for l in 0:lmax as collectPl(x; lmax)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectPl(0.5, lmax = 3)\n4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:\n  1.0\n  0.5\n -0.125\n -0.4375","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate all the associated Legendre Polynomials for coefficient m as collectPlm(x; lmax, m):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectPlm(0.5, lmax = 5, m = 3)\n6-element OffsetArray(::Vector{Float64}, 0:5) with eltype Float64 with indices 0:5:\n   0.0\n   0.0\n   0.0\n  -9.742785792574933\n -34.099750274012266\n -42.62468784251533","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate all the nth derivatives as collectdnPl(x; lmax, n):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectdnPl(0.5, lmax = 5, n = 3)\n6-element OffsetArray(::Vector{Float64}, 0:5) with eltype Float64 with indices 0:5:\n  0.0\n  0.0\n  0.0\n 15.0\n 52.5\n 65.625","category":"page"},{"location":"#Increase-precision","page":"Introduction","title":"Increase precision","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The precision of the result may be changed by using arbitrary-precision types such as BigFloat. For example, using Float64 arguments we obtain","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> Pl(1/3, 5)\n0.33333333333333337","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"whereas using BigFloat, we obtain","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> Pl(big(1)/3, 5)\n0.3333333333333333333333333333333333333333333333333333333333333333333333333333305","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The precision of the latter may be altered using setprecision, as","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> setprecision(300) do\n       Pl(big(1)/3, 5)\n       end\n0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333317","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This is particularly important to avoid overflow while computing high-order derivatives. For example:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> dnPl(0.5, 300, 200) # Float64\nNaN\n\njulia> dnPl(big(1)/2, 300, 200) # BigFloat\n1.738632750542319394663553898425873258768856732308227932150592526951212145232716e+499","category":"page"},{"location":"#Reference","page":"Introduction","title":"Reference","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Modules = [LegendrePolynomials]","category":"page"},{"location":"#LegendrePolynomials.LegendrePolynomialIterator","page":"Introduction","title":"LegendrePolynomials.LegendrePolynomialIterator","text":"LegendrePolynomialIterator(x, [lmax::Integer])\n\nReturn an iterator that generates the values of the Legendre polynomials P_ell(x) for the given x. If lmax is specified then only the values of P_ell(x) from 0 to lmax are returned.\n\nwarn: Warn\nLegendrePolynomialIterator will not be a part of the public API from the next minor release.\n\nExamples\n\njulia> import LegendrePolynomials: LegendrePolynomialIterator\n\njulia> iter = LegendrePolynomialIterator(0.5, 4);\n\njulia> collect(iter)\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  1.0\n  0.5\n -0.125\n -0.4375\n -0.2890625\n\njulia> iter = LegendrePolynomialIterator(0.5);\n\njulia> collect(Iterators.take(iter, 5)) # evaluete 5 elements (l = 0:4)\n5-element Vector{Float64}:\n  1.0\n  0.5\n -0.125\n -0.4375\n -0.2890625\n\njulia> collect(Iterators.take(Iterators.drop(iter, 100), 5)) # evaluate Pl for l = 100:104\n5-element Vector{Float64}:\n -0.0605180259618612\n  0.02196749072249231\n  0.08178451892628381\n  0.05963329258495025\n -0.021651535258316177\n\n\n\n\n\n","category":"type"},{"location":"#LegendrePolynomials.Pl-Tuple{Any, Integer}","page":"Introduction","title":"LegendrePolynomials.Pl","text":"Pl(x, l::Integer)\n\nCompute the Legendre Polynomial P_ell(x) for the argument x and the degree l.\n\nExamples\n\njulia> Pl(1, 2)\n1.0\n\njulia> Pl(0.5, 20)\n-0.04835838106737356\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.Plm","page":"Introduction","title":"LegendrePolynomials.Plm","text":"Plm(x, l::Integer, m::Integer, [cache::AbstractVector])\n\nCompute the associatedLegendre polynomial P_ell^m(x). Optionally a pre-allocated vector cache may be provided, which must have a minimum length of l - m + 1 and may be overwritten during the computation.\n\nThe polynomials are defined to include the Condon-Shortley phase (-1)^m.\n\nThe coefficient m must be non-negative. For m == 0 this function just returns Legendre polynomials.\n\nExamples\n\njulia> Plm(0.5, 3, 2)\n5.625\n\njulia> Plm(0.5, 4, 0) == Pl(0.5, 4)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"#LegendrePolynomials.collectPl!-Tuple{AbstractVector{T} where T, Any}","page":"Introduction","title":"LegendrePolynomials.collectPl!","text":"collectPl!(v::AbstractVector, x; [lmax::Integer = length(v) - 1])\n\nCompute the Legendre Polynomials P_ell(x) for the argument x and all degrees l in 0:lmax, and store them in v.\n\nAt output v[firstindex(v) + l] == Pl(x,l).\n\nExamples\n\njulia> v = zeros(4);\n\njulia> collectPl!(v, 0.5)\n4-element Vector{Float64}:\n  1.0\n  0.5\n -0.125\n -0.4375\n\njulia> v = zeros(0:4);\n\njulia> collectPl!(v, 0.5, lmax = 3) # only l from 0 to 3 are populated\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  1.0\n  0.5\n -0.125\n -0.4375\n  0.0\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPl-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectPl","text":"collectPl(x; lmax::Integer)\n\nCompute the Legendre Polynomial P_ell(x) for the argument x and all degrees l in 0:lmax. Return P with indices 0:lmax, where P[l] == Pl(x,l)\n\nExamples\n\njulia> collectPl(0.5, lmax = 4)\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  1.0\n  0.5\n -0.125\n -0.4375\n -0.2890625\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPlm!-Tuple{Any, Any}","page":"Introduction","title":"LegendrePolynomials.collectPlm!","text":"collectPlm!(v::AbstractVector, x; lmax::Integer, m::Integer)\n\nCompute the associated Legendre Polynomial P_ell^m(x) for the argument x and all degrees l = 0:lmax, and store the result in v. The polynomials for l < m are defined to be zero.\n\nThe polynomials are defined to include the Condon-Shortley phase (-1)^m.\n\nThe coefficient m must be greater than or equal to zero.\n\nAt output, v[l + firstindex(v)] == Plm(x, l, m) for l = 0:lmax.\n\nExamples\n\njulia> v = zeros(4);\n\njulia> collectPlm!(v, 0.5, lmax = 3, m = 2)\n4-element Vector{Float64}:\n 0.0\n 0.0\n 2.25\n 5.625\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPlm-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectPlm","text":"collectPlm(x; lmax::Integer, m::Integer)\n\nCompute the associated Legendre Polynomial P_ell^m(x) for the argument x and all degrees l = 0:lmax. The polynomials for l < m are defined to be zero.\n\nThe polynomials are defined to include the Condon-Shortley phase (-1)^m.\n\nThe coefficient m must be greater than or equal to zero.\n\nReturns v with indices 0:lmax, where v[l] == Plm(x, l, m).\n\nExamples\n\njulia> collectPlm(0.5, lmax = 3, m = 2)\n4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:\n 0.0\n 0.0\n 2.25\n 5.625\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectdnPl!-Tuple{Any, Any}","page":"Introduction","title":"LegendrePolynomials.collectdnPl!","text":"collectdnPl!(v::AbstractVector, x; lmax::Integer, n::Integer)\n\nCompute the n-th derivative of a Legendre Polynomial P_ell(x) for the argument x and all degrees l = 0:lmax, and store the result in v.\n\nThe order of the derivative n must be greater than or equal to zero.\n\nAt output, v[l + firstindex(v)] == dnPl(x, l, n) for l = 0:lmax.\n\nExamples\n\njulia> v = zeros(4);\n\njulia> collectdnPl!(v, 0.5, lmax = 3, n = 2)\n4-element Vector{Float64}:\n 0.0\n 0.0\n 3.0\n 7.5\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectdnPl-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectdnPl","text":"collectdnPl(x; lmax::Integer, n::Integer)\n\nCompute the n-th derivative of a Legendre Polynomial P_ell(x) for the argument x and all degrees l = 0:lmax.\n\nThe order of the derivative n must be greater than or equal to zero.\n\nReturns v with indices 0:lmax, where v[l] == dnPl(x, l, n).\n\nExamples\n\njulia> collectdnPl(0.5, lmax = 3, n = 2)\n4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:\n 0.0\n 0.0\n 3.0\n 7.5\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.dnPl","page":"Introduction","title":"LegendrePolynomials.dnPl","text":"dnPl(x, l::Integer, n::Integer, [cache::AbstractVector])\n\nCompute the n-th derivative d^n P_ell(x)dx^n of the Legendre polynomial P_ell(x). Optionally a pre-allocated vector cache may be provided, which must have a minimum length of l - n + 1 and may be overwritten during the computation.\n\nThe order of the derivative n must be non-negative. For n == 0 this function just returns Legendre polynomials.\n\nExamples\n\njulia> dnPl(0.5, 3, 2) # second derivative of P3(x) at x = 0.5\n7.5\n\njulia> dnPl(0.5, 4, 0) == Pl(0.5, 4) # zero-th order derivative == Pl(x)\ntrue\n\n\n\n\n\n","category":"function"}]
}
