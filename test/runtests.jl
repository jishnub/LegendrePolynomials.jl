using Test
using LegendrePolynomials
using OffsetArrays
using HyperDualNumbers
using Aqua

@testset "Project quality assurance" begin
    Aqua.test_all(LegendrePolynomials, ambiguities = false)
end

tohyper(x) = Hyper(x, one(x), one(x), zero(x))

@testset "Pl and collectPl" begin
	x = 2rand() - 1
    lmax = 5
    lmin = max(0, lmax - 1)

    @testset "standard norm" begin
        P = @inferred collectPl(x, lmax = lmax)
        T = Union{OffsetVector{Float64, Vector{Float64}}, OffsetVector{BigFloat, Vector{BigFloat}}}
        P2 = @inferred T collectPlm(x, lmax = lmax, m = 0)
        P3 = @inferred collectdnPl(x, lmax = lmax, n = 0)
        @test P[0] ≈ P2[0] ≈ P3[0] ≈ 1
        @test P[1] ≈ P2[1] ≈ P3[1] ≈ x
        @test P[2] ≈ P2[2] ≈ P3[2] ≈ (3x^2 - 1)/2
        @test P[3] ≈ P2[3] ≈ P3[3] ≈ (5x^3 - 3x)/2
        @test P[4] ≈ P2[4] ≈ P3[4] ≈ (35x^4 - 30x^2 + 3)/8
        @test P[5] ≈ P2[5] ≈ P3[5] ≈ (63x^5 - 70x^3 + 15x)/8

        Q = @inferred collectPl(x, lmax = lmax, lmin = lmin, norm = Val(:standard))
        Q2 = @inferred T collectPlm(x, lmax = lmax, lmin = lmin, m = 0, norm = Val(:standard))
        @test Q ≈ Q2
        @test @views Q[lmin:lmax] == P[lmin:lmax]
    end

    @testset "normalized" begin
        P = @inferred collectPl(x, lmax = lmax, norm = Val(:normalized))
        Q = @inferred collectPlm(x, lmax = lmax, m = 0, norm = Val(:normalized))
        @test P ≈ Q
        @test P[0] ≈ √(1/2)

        P2 = @inferred collectPl(x, lmax = lmax, lmin = lmin, norm = Val(:normalized))
        Q2 = @inferred collectPlm(x, lmax = lmax, lmin = lmin, m = 0, norm = Val(:normalized))
        @test P2 ≈ Q2
        @test @views P2[lmin:lmax] == P[lmin:lmax]
    end

    @testset "x = 0" begin
        for l = 1:2:101
            P1 = @inferred Pl(0, l)
            P2 = @inferred Union{BigFloat, Float64} Plm(0, l, 0)
            P3 = @inferred dnPl(0, l, 0)
            @test P1 == P2 == P3 == 0
            P1 = @inferred Pl(0, l, norm = Val(:normalized))
            P2 = @inferred Plm(0, l, 0, norm = Val(:normalized))
            @test P1 == P2 == 0
        end
    end

    @testset "x = ±1" begin
        lmax = 10

        P = collectPl(1, lmax = lmax)
        @test all(==(1), P)
        P2 = collectPlm(1, lmax = lmax, m = 0)
        @test all(≈(1), P2)
        P3 = collectdnPl(1, lmax = lmax, n = 0)
        @test all(==(1), P3)

        P = collectPl(-1, lmax = lmax)
        P2 = collectPlm(-1, lmax = lmax, m = 0)
        P3 = collectdnPl(-1, lmax = lmax, n = 0)
        for l in axes(P, 1)
            @test P[l] ≈ P2[l] ≈ P3[l] ≈ (-1)^l
        end
    end

    @test_throws DomainError collectPl(-2, lmax = lmax)
    @test_throws ArgumentError collectPl(0, lmax = -1)
end

@testset "Plm" begin
    @testset "double factorial" begin
        @test Plm(0.5, 2, 2) ≈ 9/4
        @test Plm(0.5, 3, 3) ≈ -45/8*sqrt(3)
        @test Plm(0.5, 3, 3) ≈ -45/8*sqrt(3)
        @test Plm(0.5, 10, 0) ≈ -49343/262144
        @test Plm(0.5, 10, 6) ≈ -674999325/8192
        @test Plm(0.5, 10, 10) ≈ 159099165225/1024
    end
    @test_throws ArgumentError Plm(0.5, 3, 10)
    @test_throws ArgumentError Plm(0.5, -1, 3)
    @test_throws ArgumentError Plm(0.5, 1, -3)
    @test_throws DomainError Plm(2.5, 1, 3)
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

@testset "collectPlm!" begin
    v = zeros(10);
    @test_throws DomainError collectPlm!(v, -2, lmax = 3, m = 2)
    @test_throws ArgumentError collectPlm!(v, 0, lmax = -1, m = 3)
    @test_throws ArgumentError collectPlm!(v, 0.5, lmax = 1, m = -3)
end

@testset "collectPlm" begin
    @test_throws DomainError collectPlm(-2, lmax = 3, m = 2)
    @test_throws ArgumentError collectPlm(0, lmax = -1, m = 3)
    @test_throws ArgumentError collectPlm(0.5, lmax = 1, m = -3)
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
    P2 = collectPlm(x, lmax = lmax, m = 1) / (- sqrt(1 - x^2))
    P3 = collectdnPl(x, lmax = lmax, n = 1)
    @test P[0] == P3[0] == 0
    @test P[1] ≈ P2[1] ≈ P3[1] ≈ 1
    @test P[2] ≈ P2[2] ≈ P3[2] ≈ 3x
    @test P[3] ≈ P2[3] ≈ P3[3] ≈ (-3 + 15x^2)/2
    @test P[4] ≈ P2[4] ≈ P3[4] ≈ 1/8*(-60x + 140x^3)
    @test P[5] ≈ P2[5] ≈ P3[5] ≈ 1/8*(15 - 210x^2 + 315x^4)

    @testset "x = ±1" begin
        lmax = 10
        x = 1
        P = collectPl(tohyper(x), lmax = lmax)
        dP1h = eps1.(P)
        dP1 = collectdnPl(x, lmax = lmax, n = 1)

        for l in 1:lmax
            @test dP1h[l] == dP1[l] == l*(l+1)/2
        end

        x = -1
        P = collectPl(tohyper(x), lmax = lmax)
        dPm1h = eps1.(P)
        dPm1 = collectdnPl(x, lmax = lmax, n = 1)

        for l in 1:lmax
            @test dPm1h[l] == dPm1[l] == (-1)^(l+1) * l*(l+1)/2
        end
    end
end

@testset "d^2Pl/dx^2" begin
    x = 2rand() - 1
    xh = tohyper(x)
    lmax = 5
    P = eps1eps2.(collectPl(xh, lmax = lmax))
    P2 = collectPlm(x, lmax = lmax, m = 2) / (1 - x^2)
    P3 = collectdnPl(x, lmax = lmax, n = 2)

    @test P[0] == P3[0] == 0
    @test P[1] == P3[1] == 0
    @test P[2] ≈ P2[2] ≈ P3[2] ≈ 3
    @test P[3] ≈ P2[3] ≈ P3[3] ≈ 15x
    @test P[4] ≈ P2[4] ≈ P3[4] ≈ 1/8*(-60 + 420x^2)
    @test P[5] ≈ P2[5] ≈ P3[5] ≈ 1/8*(-420x + 1260x^3)
end

@testset "Askey-Gasper" begin
    x = 0.5
    for lmax=1:100:1_000
        @test sum(collectPl(x, lmax = lmax)) >= 0
        @test sum(collectPlm(x, lmax = lmax, m = 0)) >= 0
        @test sum(collectdnPl(x, lmax = lmax, n = 0)) >= 0
    end
end
