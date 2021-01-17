using Test
using LegendrePolynomials
using OffsetArrays
using HyperDualNumbers
using Aqua

@testset "Project quality assurance" begin
    Aqua.test_all(LegendrePolynomials)
end

import LegendrePolynomials: LegendrePolynomialIterator

tohyper(x) = Hyper(x, one(x), one(x), zero(x))

@testset "LegendrePolynomialIterator" begin
    lmax = 5
    iter = LegendrePolynomialIterator(0.5, lmax);
    @test length(iter) == lmax + 1
    @test eltype(iter) == Float64
    @test eltype(typeof(iter)) == Float64
    @test size(iter) == (lmax+1,)
    @test axes(iter) == (0:lmax,)
    @test keys(iter) == 0:lmax

    iter2 = copy(iter)
    @test typeof(iter2) == typeof(iter)
    @test iter2.x == iter.x
    @test iter2.lmax == iter.lmax

    iter = LegendrePolynomialIterator(0.5);
    @test Base.IteratorSize(typeof(iter)) == Base.IsInfinite()

    @test_throws DomainError LegendrePolynomialIterator(-1.1, lmax)
    @test_throws DomainError LegendrePolynomialIterator(-1.1)
end

@testset "Pl" begin
	x = 2rand() - 1 
    lmax = 5
    P = collectPl(x, lmax = lmax)
    P2 = collectdnPl(x, lmax = lmax, n = 0)
    @test P[0] == P2[0] == 1
    @test P[1] == P2[1] == x
    @test P[2] ≈ P2[2] ≈ (3x^2 - 1)/2
    @test P[3] ≈ P2[3] ≈ (5x^3 - 3x)/2
    @test P[4] ≈ P2[4] ≈ (35x^4 - 30x^2 + 3)/8
    @test P[5] ≈ P2[5] ≈ (63x^5 - 70x^3 + 15x)/8

    @testset "x = 0" begin
        for l = 1:2:101
            @test Pl(0, l) == dnPl(0, l, 0) == 0
        end
    end

    @testset "x = ±1" begin
        lmax = 10

        P = collectPl(1, lmax = lmax)
        @test all(==(1), P)
        P2 = collectdnPl(1, lmax = lmax, n = 0)
        @test all(==(1), P2)

        P = collectPl(-1, lmax = lmax)
        P2 = collectdnPl(-1, lmax = lmax, n = 0)
        for l in axes(P, 1)
            @test P[l] == P2[l] == (-1)^l
        end
    end

    @test_throws DomainError collectPl(-2, lmax = lmax)
    @test_throws ArgumentError collectPl(0, lmax = -1)
end

@testset "dnPl" begin
    @test dnPl(0.5, 3, 10) == 0
    @test dnPl(0.5, 3, 10, zeros(100)) == 0
    @testset "double factorial" begin
        @test dnPl(0.5, 2, 2) == 3
        @test dnPl(0.5, 3, 3) == 15
        @test dnPl(0.5, 3, 3, zeros(1)) == 15
    end
    @test_throws ArgumentError dnPl(0.5, 10, 3, zeros(1))
    @test_throws ArgumentError dnPl(0.5, -1, 3)
    @test_throws ArgumentError dnPl(0.5, 1, -3)
    @test_throws DomainError dnPl(2.5, 1, 3)
end

@testset "collectPl!" begin
    v = zeros(10);
    x = 0.6

    collectPl!(v, x)
    for l in 0:length(v) - 1
        @test v[firstindex(v) + l] == Pl(x, l)
    end

    v .= 0
    lmax = 6
    collectPl!(v, x, lmax = lmax)
    for l in 0:lmax
        @test v[firstindex(v) + l] == Pl(x, l)
    end

    @test_throws DomainError collectPl!(v, -2, lmax = lmax)
    @test_throws ArgumentError collectPl!(v, 0, lmax = -1)
end

@testset "collectdnPl!" begin
    v = zeros(10);
    @test_throws DomainError collectdnPl!(v, -2, lmax = 3, n = 2)
    @test_throws ArgumentError collectdnPl!(v, 0, lmax = -1, n = 3)
    @test_throws ArgumentError collectdnPl!(v, 0.5, lmax = 1, n = -3)
    @test_throws ArgumentError collectdnPl!(zeros(1), 0, lmax = 3, n = 3)
end

@testset "collectdnPl" begin
    @test_throws DomainError collectdnPl(-2, lmax = 3, n = 2)
    @test_throws ArgumentError collectdnPl(0, lmax = -1, n = 3)
    @test_throws ArgumentError collectdnPl(0.5, lmax = 1, n = -3)
end

@testset "dPl/dx" begin
    x = 2rand() - 1 
    xh = tohyper(x)
    lmax = 5
    P = eps1.(collectPl(xh, lmax = lmax))
    P2 = collectdnPl(x, lmax = lmax, n = 1)
    @test P[0] == P2[0] == 0
    @test P[1] == P2[1] == 1
    @test P[2] ≈ P2[2] ≈ 3x
    @test P[3] ≈ P2[3] ≈ (-3 + 15x^2)/2
    @test P[4] ≈ P2[4] ≈ 1/8*(-60x + 140x^3)
    @test P[5] ≈ P2[5] ≈ 1/8*(15 - 210x^2 + 315x^4)

    @testset "x = ±1" begin
        lmax = 10
        x = 1
        P = collectPl(tohyper(x), lmax = lmax)
        dP1h = eps1.(P)
        dP1 = collectdnPl(x, lmax = lmax, n = 1)

        for l in axes(dP1,1)
            @test dP1h[l] == dP1[l] == l*(l+1)/2
        end

        x = -1
        P = collectPl(tohyper(x), lmax = lmax)
        dPm1h = eps1.(P)
        dPm1 = collectdnPl(x, lmax = lmax, n = 1)

        for l in axes(dPm1,1)
            @test dPm1h[l] == dPm1[l] == (-1)^(l+1) * l*(l+1)/2
        end
    end
end

@testset "d^2Pl/dx^2" begin
    x = 2rand() - 1 
    xh = tohyper(x)
    lmax = 5
    P = eps1eps2.(collectPl(xh, lmax = lmax))
    P2 = collectdnPl(x, lmax = lmax, n = 2)

    @test P[0] == P2[0] == 0
    @test P[1] == P2[1] == 0

    @test P[2] ≈ P2[2] ≈ 3
    @test P[3] ≈ P2[3] ≈ 15x
    @test P[4] ≈ P2[4] ≈ 1/8*(-60 + 420x^2)
    @test P[5] ≈ P2[5] ≈ 1/8*(-420x + 1260x^3)
end

@testset "Askey-Gasper" begin
    x = 0.5
    for lmax=1:100:1_000
        @test sum(collectPl(x, lmax = lmax)) >= 0
        @test sum(collectdnPl(x, lmax = lmax, n = 0)) >= 0
    end
end