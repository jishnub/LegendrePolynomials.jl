using Test
using LegendrePolynomials
using OffsetArrays
using HyperDualNumbers

tohyper(x) = Hyper(x, one(x), one(x), zero(x))

@testset "Pl" begin
	x = 2rand() - 1 
    P = collectPl(x, lmax = 5)
    @test P[0] ≈ 1
    @test P[1] ≈ x
    @test P[2] ≈ (3x^2 - 1)/2
    @test P[3] ≈ (5x^3 - 3x)/2
    @test P[4] ≈ (35x^4 - 30x^2 + 3)/8
    @test P[5] ≈ (63x^5 - 70x^3 + 15x)/8

    @testset "x = 0" begin
        for l = 1:2:101
            @test Pl(0, l) == 0
        end
    end

    @testset "x = ±1" begin
        lmax = 10

        P = collectPl(1, lmax = lmax)
        @test all(==(1), P)

        P = collectPl(-1, lmax = lmax)
        for l in axes(P, 1)
            @test P[l] == (-1)^l
        end
    end
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
end

@testset "dPl/dx" begin
    x = 2rand() - 1 
    xh = tohyper(x)
    P = eps1.(collectPl(xh, lmax=5))
    @test P[0] ≈ 0
    @test P[1] ≈ 1
    @test P[2] ≈ 3x
    @test P[3] ≈ (-3 + 15x^2)/2
    @test P[4] ≈ 1/8*(-60x + 140x^3)
    @test P[5] ≈ 1/8*(15 - 210x^2 + 315x^4)

    @testset "x = ±1" begin
        lmax = 10
        
        P = collectPl(tohyper(1), lmax = lmax)
        dP1 = eps1.(P)
        for l in axes(dP1,1)
            @test dP1[l] == l*(l+1)/2
        end

        P = collectPl(tohyper(-1), lmax = lmax)
        dPm1 = eps1.(P)
        for l in axes(dPm1,1)
            @test dPm1[l] == (-1)^(l+1) * l*(l+1)/2
        end
    end
end

@testset "d^2Pl/dx^2" begin
    x = 2rand() - 1 
    xh = tohyper(x)
    P = eps1eps2.(collectPl(xh, lmax=5))
    @test P[0] ≈ 0
    @test P[1] ≈ 0
    @test P[2] ≈ 3
    @test P[3] ≈ 15x
    @test P[4] ≈ 1/8*(-60 + 420x^2)
    @test P[5] ≈ 1/8*(-420x + 1260x^3)
end

@testset "Askey-Gasper" begin
    x = 0.5
    for lmax=1:100:1_000
        @test sum(collectPl(x, lmax = lmax)) >= 0
    end
end