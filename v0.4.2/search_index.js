var documenterSearchIndex = {"docs":
[{"location":"derivatives/#Derivatives-of-Legendre-Polynomials","page":"Derivatives","title":"Derivatives of Legendre Polynomials","text":"","category":"section"},{"location":"derivatives/#Analytical-recursive-approach","page":"Derivatives","title":"Analytical recursive approach","text":"","category":"section"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Bonnet's recursion formula","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"P_ell(x) = left((2ell-1) x P_ell-1(x) - (ell-1)P_ell - 2(x)right)ell","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"may be differentiated an arbitrary number of times analytically to obtain recursion relations for higher derivatives:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"fracd^n P_ell(x)dx^n = frac(2ell-1)ell left(x fracd^n P_ell-1(x)dx^n +\nn fracd^(n-1) P_ell-1(x)dx^(n-1) right) - frac(ell-1)ell fracd^n P_ell-2(x)dx^n","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"This provides a simultaneous recursion relation in ell as well as n, solving which we may obtain derivatives up to any order. This is the approach used in this package to compute the derivatives of Legendre polynomials.","category":"page"},{"location":"derivatives/#Automatic-diferentiation","page":"Derivatives","title":"Automatic diferentiation","text":"","category":"section"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Julia automatic differentiation framework may be used to compute the derivatives of Legendre polynomials alongside their values. Since the defintions of the polynomials are completely general, they may be called with dual or hyperdual numbers as arguments to evaluate derivarives in one go. We demonstrate one example of this using the package HyperDualNumbers.jl v4:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"DocTestSetup = quote\n\tusing LegendrePolynomials\n  using HyperDualNumbers\nend","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> x = 0.5;\n\njulia> Pl(x, 3)\n-0.4375\n\njulia> using HyperDualNumbers\n\njulia> xh = Hyper(x, one(x), one(x), zero(x));\n\njulia> p = Pl(xh, 3)\n-0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The Legendre polynomial P_ell(x) may be obtained using","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> realpart(p)\n-0.4375","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The first derivative dP_ell(x)dx may be obtained as","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁part(p)\n0.375","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"The second derivative d^2P_ell(x)dx^2 may be obtained using","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁ε₂part(p)\n7.5","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Something similar may also be evaluated for all l iteratively. For example, the function collectPl evaluates Legendre polynomials for the HyperDualNumber argument for a range of l:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> collectPl(xh, lmax=4)\n5-element OffsetArray(::Vector{Hyper{Float64}}, 0:4) with eltype Hyper{Float64} with indices 0:4:\n                1.0 + 0.0ε₁ + 0.0ε₂ + 0.0ε₁ε₂\n                0.5 + 1.0ε₁ + 1.0ε₂ + 0.0ε₁ε₂\n             -0.125 + 1.5ε₁ + 1.5ε₂ + 3.0ε₁ε₂\n        -0.4375 + 0.375ε₁ + 0.375ε₂ + 7.5ε₁ε₂\n -0.2890625 - 1.5625ε₁ - 1.5625ε₂ + 5.625ε₁ε₂","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"We may extract the first derivatives by broadcasting the function ε₁part on the array as:","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> ε₁part.(collectPl(xh, lmax=4))\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  0.0\n  1.0\n  1.5\n  0.375\n -1.5625","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Similarly the function ε₁ε₂part may be used to obtain the second derivatives.","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Several convenience functions to compute the derivatives of Legendre polynomials were available in LegendrePolynomials v0.2, but have been removed in v0.3. The users are encouraged to implement convenience functions to extract the derivatives as necessary. As an exmaple, we may compute the polynomials and their first and second derivatives together as","category":"page"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"julia> using HyperDualNumbers\n\njulia> function Pl_dPl_d2Pl(x; lmax)\n           xh = Hyper(x, one(x), one(x), zero(x))\n           p = collectPl(xh, lmax = lmax)\n           realpart.(p), ε₁part.(p), ε₁ε₂part.(p)\n       end\nPl_dPl_d2Pl (generic function with 1 method)\n\njulia> Pl_dPl_d2Pl(0.5, lmax = 3)\n([1.0, 0.5, -0.125, -0.4375], [0.0, 1.0, 1.5, 0.375], [0.0, 0.0, 3.0, 7.5])","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = LegendrePolynomials\nDocTestSetup = quote\n\tusing LegendrePolynomials\nend","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Compute Legendre polynomials, Associated Legendre polynomials, and their derivatives using a 3-term recursion relation (Bonnet’s recursion formula).","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"P_ell(x) = left((2ell-1) x P_ell-1(x) - (ell-1)P_ell - 2(x)right)ell","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"By default this package evaluates the standard polynomials that satisfy P_ell(1) = 1 and P_0(x) = 1. These are normalized as","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"int_-1^1 P_m(x) P_n(x) dx = frac22n+1 delta_mn","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Optionally, normalized polynomials may be evaluated that have an L2 norm of 1.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Analogous to Legendre polynomials, one may evaluate associated Legendre polynomials using a 3-term recursion relation. This is evaluated by iterating over the normalized associated Legendre functions, and multiplying the norm at the final stage. Such an iteration avoids floating-point overflow.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The relation used in evaluating normalized associated Legendre polynomials barP_ell^mleft(xright) is","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"barP_ell^mleft(xright)=alpha_ell mleft(xbarP_ell-1^mleft(xright)-frac1alpha_ell-1mbarP_ell-2^mleft(xright)right)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"where","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"alpha_ell m=sqrtfracleft(2ell+1right)left(2ell-1right)left(ell-mright)left(ell+mright)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The starting value for a specific m is obtained by setting ell = m, in which case we use the explicit relation","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"barP_m^mleft(xright)=left(-1right)^mleft(2m-1right)sqrtfracleft(2m+1right)2left(2mright)left(1-x^2right)^left(m2right)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The polynomials, thus defined, include the Condon-Shortley phase (-1)^m, and may be evaluated using the standard normalization as","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"P_ell^mleft(xright)=sqrtfrac2left(ell+mright)left(2ell+1right)left(ell-mright)barP_ell^mleft(xright)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"There are six main functions:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pl(x,l; [norm = Val(:standard)]): this evaluates the Legendre polynomial for a given degree l at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectPl(x; lmax, [lmin = 0], [norm = Val(:standard)]): this evaluates all the polynomials for l lying in 0:lmax at the argument x. As before the argument needs to lie in the domain of validity. Functionally this is equivalent to Pl.(x, lmin:lmax), except collectPl evaluates the result in one pass, and is therefore faster. There is also the in-place version collectPl! that uses a pre-allocated array.\nPlm(x, l, m; [norm = Val(:standard)]): this evaluates the associated Legendre polynomial P_ell^m(x) at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectPlm(x; m, lmax, [lmin = abs(m)], [norm = Val(:standard)]): this evaluates the associated Legendre polynomials with coefficient m for l = lmin:lmax. There is also an in-place version collectPlm! that uses a pre-allocated array.\ndnPl(x, l, n): this evaluates the n-th derivative of the Legendre polynomial P_ell(x) at the argument x. The argument needs to satisfy -1 <= x <= 1.\ncollectdnPl(x; n, lmax): this evaluates the n-th derivative of all the Legendre polynomials for l = 0:lmax. There is also an in-place version collectdnPl! that uses a pre-allocated array.","category":"page"},{"location":"#Quick-Start","page":"Introduction","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the Legendre polynomial for one l at an argumentx as Pl(x, l):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> p = Pl(0.5, 3)\n-0.4375\n\njulia> p ≈ -7/16 # analytical value\ntrue","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the associated Legendre Polynomial one l,m pair as Plm(x, l, m):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> p = Plm(0.5, 3, 2)\n5.624999999999997\n\njulia> p ≈ 45/8 # analytical value\ntrue","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the nth derivative of Pl for one l as dnPl(x, l, n):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> dnPl(0.5, 3, 2)\n7.5","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the Legendre polynomials for l in lmin:lmax as collectPl(x; lmin, lmax). By default lmin is chosen to be 0, and may be omitted.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectPl(0.5, lmax = 3)\n4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:\n  1.0\n  0.5\n -0.125\n -0.4375","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate the associated Legendre Polynomials for order m and degree l in lmin:lmax as collectPlm(x; m, lmin, lmax). By default lmin is chosen to be abs(m), and may be omitted.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectPlm(0.5, lmax = 5, m = 3)\n3-element OffsetArray(::Vector{Float64}, 3:5) with eltype Float64 with indices 3:5:\n  -9.742785792574933\n -34.09975027401223\n -42.62468784251535","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Evaluate all the nth derivatives of Pl as collectdnPl(x; lmax, n):","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> collectdnPl(0.5, lmax = 5, n = 3)\n6-element OffsetArray(::Vector{Float64}, 0:5) with eltype Float64 with indices 0:5:\n  0.0\n  0.0\n  0.0\n 15.0\n 52.5\n 65.625","category":"page"},{"location":"#Reference","page":"Introduction","title":"Reference","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Modules = [LegendrePolynomials]","category":"page"},{"location":"#LegendrePolynomials.Pl-Tuple{Any, Integer}","page":"Introduction","title":"LegendrePolynomials.Pl","text":"Pl(x, l::Integer; [norm = Val(:standard)])\n\nCompute the Legendre Polynomial P_ell(x) for the argument x and the degree l.\n\nThe default norm is chosen to be Val(:standard), in which case the polynomials satisfy Pl(1, l) == 1 for all l. Optionally, the norm may be set to Val(:normalized) to obtain normalized Legendre polynomials. These have an L2 norm of 1.\n\nExamples\n\njulia> Pl(1, 2)\n1.0\n\njulia> Pl(0.5, 4) ≈ -37/128 # analytically obtained value\ntrue\n\njulia> Pl(0.5, 20, norm = Val(:normalized))\n-0.21895188261094017\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.Plm-Tuple{Any, Integer, Integer}","page":"Introduction","title":"LegendrePolynomials.Plm","text":"Plm(x, l::Integer, m::Integer; [kwargs...])\n\nCompute the associated Legendre polynomial P_ell^m(x) for degree l and order m lying in -l:l.\n\nFor m == 0 this function returns Legendre polynomials.\n\nOptional keyword arguments\n\nnorm: The normalization used in the function. Possible   options are Val(:standard) (default) and Val(:normalized).   The functions obtained with norm = Val(:normalized) have an L2 norm of 1.\ncsphase::Bool: The Condon-shortley phase (-1)^m, which is included by default.\n\nExamples\n\njulia> Plm(0.5, 3, 2) ≈ 45/8 # analytically obtained value\ntrue\n\njulia> Plm(0.5, 4, 0) == Pl(0.5, 4)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPl!-Tuple{AbstractVector{T} where T, Any}","page":"Introduction","title":"LegendrePolynomials.collectPl!","text":"collectPl!(v::AbstractVector, x; [lmin::Integer = 0], [lmax::Integer = length(v) - 1 + lmin], [kwargs...])\n\nCompute the Legendre Polynomials P_ell(x) for the argument x and degrees l = lmin:lmax, and store them in v.\n\nAt output v[firstindex(v) + l] == Pl(x, l; kwargs...).\n\nOptional keyword arguments\n\nnorm: The normalization used in the function. Possible   options are Val(:standard) (default) and Val(:normalized).   The functions obtained with norm = Val(:normalized) have an L2 norm of 1.\n\nExamples\n\njulia> v = zeros(4);\n\njulia> collectPl!(v, 0.5);\n\njulia> all(zip(0:3, v)) do (l, vl); vl ≈ Pl(0.5, l); end\ntrue\n\njulia> collectPl!(v, 0.5, norm = Val(:normalized));\n\njulia> all(zip(0:3, v)) do (l,vl); vl ≈ Pl(0.5, l, norm = Val(:normalized)); end\ntrue\n\njulia> v = zeros(0:4);\n\njulia> collectPl!(v, 0.5, lmax = 3) # only l from 0 to 3 are populated\n5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:\n  1.0\n  0.5\n -0.125\n -0.4375\n  0.0\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPl-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectPl","text":"collectPl(x; lmax::Integer, [lmin::Integer = 0], [kwargs...])\n\nCompute the Legendre Polynomial P_ell(x) for the argument x and degrees l = lmin:lmax. Return P with indices lmin:lmax, where P[l] == Pl(x, l; kwargs...)\n\nOptional keyword arguments\n\nnorm: The normalization used in the function. Possible   options are Val(:standard) (default) and Val(:normalized).   The functions obtained with norm = Val(:normalized) have an L2 norm of 1.\n\nExamples\n\njulia> v = collectPl(0.5, lmax = 4);\n\njulia> all(l -> v[l] ≈ Pl(0.5, l), 0:4)\ntrue\n\njulia> v = collectPl(0.5, lmax = 4, norm = Val(:normalized));\n\njulia> all(l -> v[l] ≈ Pl(0.5, l, norm = Val(:normalized)), 0:4)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPlm!-Tuple{Any, Any}","page":"Introduction","title":"LegendrePolynomials.collectPlm!","text":"collectPlm!(v::AbstractVector, x; m::Integer, [lmin::Integer = abs(m)], [lmax::Integer = length(v) - 1 + lmin], [kwargs...])\n\nCompute the associated Legendre polynomials P_ell^m(x) for the argument x, degrees l = lmin:lmax and order m lying in -l:l, and store the result in v.\n\nAt output, v[l + firstindex(v)] == Plm(x, l, m; kwargs...) for l = lmin:lmax.\n\nOptional keyword arguments\n\nnorm: The normalization used in the function. Possible   options are Val(:standard) (default) and Val(:normalized).   The functions obtained with norm = Val(:normalized) have an L2 norm of 1.\ncsphase::Bool: The Condon-shortley phase (-1)^m, which is included by default.\n\nExamples\n\njulia> v = zeros(2);\n\njulia> collectPlm!(v, 0.5, lmax = 3, m = 2);\n\njulia> all(zip(2:3, v)) do (l, vl); vl ≈ Plm(0.5, l, 2); end\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectPlm-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectPlm","text":"collectPlm(x; m::Integer, lmax::Integer, [lmin::Integer = abs(m)], [kwargs...])\n\nCompute the associated Legendre polynomials P_ell^m(x) for the argument x, degrees l = lmin:lmax and order m lying in -l:l.\n\nReturns v with indices lmin:lmax, with v[l] == Plm(x, l, m; kwargs...).\n\nOptional keyword arguments\n\nnorm: The normalization used in the function. Possible   options are Val(:standard) (default) and Val(:normalized).   The functions obtained with norm = Val(:normalized) have an L2 norm of 1.\ncsphase::Bool: The Condon-shortley phase (-1)^m, which is included by default.\nTnorm: Type of the normalization factor, which is dynamicall inferred by default,   and is used in allocating an appropriate container.   In case it is known that the norm may be expressed as a certain type without overflow   (eg. Float64), this may be provided. Care must be taken to avoid overflows if Tnorm   is passed as an argument.\n\nExamples\n\njulia> v = collectPlm(0.5, lmax = 3, m = 2);\n\njulia> all(l -> v[l] ≈ Plm(0.5, l, 2), 2:3)\ntrue\n\njulia> v = collectPlm(0.5, lmax = 3, m = 2, norm = Val(:normalized));\n\njulia> all(l -> v[l] ≈ Plm(0.5, l, 2, norm = Val(:normalized)), 2:3)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectdnPl!-Tuple{Any, Any}","page":"Introduction","title":"LegendrePolynomials.collectdnPl!","text":"collectdnPl!(v::AbstractVector, x; lmax::Integer, n::Integer)\n\nCompute the n-th derivative of a Legendre Polynomial P_ell(x) for the argument x and all degrees l = 0:lmax, and store the result in v.\n\nThe order of the derivative n must be greater than or equal to zero.\n\nAt output, v[l + firstindex(v)] == dnPl(x, l, n) for l = 0:lmax.\n\nExamples\n\njulia> v = zeros(4);\n\njulia> collectdnPl!(v, 0.5, lmax = 3, n = 2)\n4-element Vector{Float64}:\n 0.0\n 0.0\n 3.0\n 7.5\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.collectdnPl-Tuple{Any}","page":"Introduction","title":"LegendrePolynomials.collectdnPl","text":"collectdnPl(x; lmax::Integer, n::Integer)\n\nCompute the n-th derivative of a Legendre Polynomial P_ell(x) for the argument x and all degrees l = 0:lmax.\n\nThe order of the derivative n must be greater than or equal to zero.\n\nReturns v with indices 0:lmax, where v[l] == dnPl(x, l, n).\n\nExamples\n\njulia> collectdnPl(0.5, lmax = 3, n = 2)\n4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:\n 0.0\n 0.0\n 3.0\n 7.5\n\n\n\n\n\n","category":"method"},{"location":"#LegendrePolynomials.dnPl","page":"Introduction","title":"LegendrePolynomials.dnPl","text":"dnPl(x, l::Integer, n::Integer, [cache::AbstractVector])\n\nCompute the n-th derivative d^n P_ell(x)dx^n of the Legendre polynomial P_ell(x). Optionally a pre-allocated vector cache may be provided, which must have a minimum length of l - n + 1 and may be overwritten during the computation.\n\nThe order of the derivative n must be non-negative. For n == 0 this function just returns Legendre polynomials.\n\nExamples\n\njulia> dnPl(0.5, 3, 2) # second derivative of P3(x) at x = 0.5\n7.5\n\njulia> dnPl(0.5, 4, 0) == Pl(0.5, 4) # zero-th order derivative == Pl(x)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"DocTestSetup = :(using LegendrePolynomials)","category":"page"},{"location":"tutorial/#Computing-normalized-Legendre-polynomials","page":"Tutorial","title":"Computing normalized Legendre polynomials","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default, the Legendre and associated Legendre polynomials are not normalized. One may specify the normalization through the keyword argument norm. The normalization options accepted are","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"norm = Val(:standard): standard, unnormalized polynomials. This is the default option.\nnorm = Val(:normalized): fully normalized polynomials with an L2 norm of 1","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"note: Note\nIrrespective of the norm specified, the 3-term recursion relations used are stable ones, so polynomials of high degrees may be computed without overflow.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> l = 3;\n\njulia> Pl(0.5, l)\n-0.4375\n\njulia> Pl(0.5, l, norm = Val(:normalized))\n-0.8184875533567997\n\njulia> Pl(0.5, l, norm = Val(:normalized)) ≈ √((2l+1)/2) * Pl(0.5, l)\ntrue\n\njulia> l = m = 3000;\n\njulia> Plm(0.5, l, m, norm = Val(:normalized))\n2.172276347346834e-187","category":"page"},{"location":"tutorial/#Condon-Shortley-phase","page":"Tutorial","title":"Condon-Shortley phase","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default, the associated Legendre polynomials are computed by including the Condon-Shortley phase (-1)^m. This may be disabled by passing the flag csphase = false to the various Plm functions.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"note: Note\nThe Pl functions do not accept the keyword argument csphase.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> Plm(0.5, 3, 3)\n-9.742785792574933\n\njulia> Plm(0.5, 3, 3, csphase = false)\n9.742785792574933\n\njulia> all(Plm(0.5, 3, m, csphase = false) ≈ (-1)^m * Plm(0.5, 3, m) for m in -3:3)\ntrue","category":"page"},{"location":"tutorial/#Polynomials-at-multiple-points","page":"Tutorial","title":"Polynomials at multiple points","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The collectPl(m) functions return a vector of polynomials evaluated for one theta. To evaluate the polynomials for multiple theta in one go, one may run","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> θr = range(0, pi, length=6);\n\njulia> collectPl.(cos.(θr), lmax=3)\n6-element Vector{OffsetArrays.OffsetVector{Float64, Vector{Float64}}}:\n [1.0, 1.0, 1.0, 1.0]\n [1.0, 0.8090169943749475, 0.4817627457812106, 0.1102457514062632]\n [1.0, 0.30901699437494745, -0.35676274578121053, -0.38975424859373686]\n [1.0, -0.30901699437494734, -0.35676274578121064, 0.3897542485937368]\n [1.0, -0.8090169943749473, 0.48176274578121037, -0.11024575140626285]\n [1.0, -1.0, 1.0, -1.0]","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This returns a vector of vectors, where each element corresponds to one theta. Often, one wants a Matrix where each column corresponds to one theta. We may obtain this as","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> mapreduce(hcat, θr) do θ\n           collectPl(cos(θ), lmax=3)\n       end\n4×6 Matrix{Float64}:\n 1.0  1.0        1.0        1.0        1.0        1.0\n 1.0  0.809017   0.309017  -0.309017  -0.809017  -1.0\n 1.0  0.481763  -0.356763  -0.356763   0.481763   1.0\n 1.0  0.110246  -0.389754   0.389754  -0.110246  -1.0","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"As of Julia v1.7.2 and OffsetArrays v1.10.8, this strips off the axis information, so one would need to wrap the result in an OffsetArray again to index the Matrix using the degrees ell.","category":"page"},{"location":"tutorial/#Increasing-precision","page":"Tutorial","title":"Increasing precision","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The precision of the result may be changed by using arbitrary-precision types such as BigFloat. For example, using Float64 arguments we obtain","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> Pl(1/3, 5)\n0.33333333333333337","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"whereas using BigFloat, we obtain","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> Pl(big(1)/3, 5)\n0.3333333333333333333333333333333333333333333333333333333333333333333333333333305","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This is particularly important to avoid overflow while computing high-order derivatives. For example:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> dnPl(0.5, 300, 200) # Float64\nNaN\n\njulia> dnPl(big(1)/2, 300, 200) # BigFloat\n1.738632750542319394663553898425873258768856732308227932150592526951212145232716e+499","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In general, one would need to use higher precision for both the argument x and the degrees l to obtain accurate results. For example:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> Plm(big(1)/2, big(3000), 3000)\n2.05451584347939475644802290993338963448971107938391335496027846832433343889916e+9844","category":"page"},{"location":"tutorial/#Symbolic-evaluation","page":"Tutorial","title":"Symbolic evaluation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"It's possible to symbolically evaluate Legendre or associated Legendre polynomials using Symbolics.jl. An exmaple is:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> using Symbolics\n\njulia> @variables x;\n\njulia> Pl(x, 3)\n(5//3)*x*((3//2)*(x^2) - (1//2)) - (2//3)*x\n\njulia> myf = eval(build_function(Pl(x, 3), [x]));\n\njulia> myf(0.4) == Pl(0.4, 3)\ntrue","category":"page"}]
}
