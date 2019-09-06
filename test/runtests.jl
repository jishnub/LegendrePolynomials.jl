using OffsetArrays,Test
using LegendrePolynomials

@testset "Pl array" begin
	x = 2rand() - 1 
    P = Pl(x,lmax=5)
    @test P[0]≈1
    @test P[1]≈x
    @test P[2]≈(3x^2-1)/2
    @test P[3]≈(5x^3-3x)/2
    @test P[4]≈(35x^4-30x^2+3)/8
    @test P[5]≈(63x^5-70x^3+15x)/8
end

@testset "dPl array" begin
	x = 2rand() - 1 
    P = dPl(x,lmax=5)
    @test P[0]≈0
    @test P[1]≈1
    @test P[2]≈3x
    @test P[3]≈(-3+15x^2)/2
    @test P[4]≈1/8*(-60x + 140x^3)
    @test P[5]≈1/8*(15 - 210x^2 + 315x^4)
end

@testset "d2Pl array" begin
	x = 2rand() - 1 
    P = d2Pl(x,lmax=5)
    @test P[0]≈0
    @test P[1]≈0
    @test P[2]≈3
    @test P[3]≈15x
    @test P[4]≈1/8*(-60 + 420x^2)
    @test P[5]≈1/8*(-420x + 1260x^3)
end

@testset "allmodes vs one" begin
    x = 2rand() - 1
    lmax = 10
    P,dP,d2P = Pl_dPl_d2Pl(x,lmax=lmax)
    @testset "Pl" begin
        for l=0:lmax
            @test P[l]==Pl(x,l)
        end
    end
    @testset "dPl" begin
        for l=0:lmax
            @test dP[l]==dPl(x,l)
        end
    end
    @testset "d2Pl" begin
        for l=0:lmax
            @test d2P[l]==d2Pl(x,l)
        end
    end 
end

@testset "boundary" begin

    lmax = 10
    P1,dP1 = Pl_dPl(1,lmax=lmax)
    Pm1,dPm1 = Pl_dPl(-1,lmax=lmax)

    @testset "Pl" begin
        for l in axes(P1,1)
            @test P1[l] == 1
        end
        for l in axes(Pm1,1)
            @test Pm1[l] == (-1)^l
        end
    end
    @testset "dPl" begin
        for l in axes(dP1,1)
            @test dP1[l] == l*(l+1)/2
        end
        for l in axes(dPm1,1)
            @test dPm1[l] == (-1)^(l+1) * l*(l+1)/2
        end
    end
end

@testset "Askey-Gasper" begin
    x=0.5
    for lmax=1:100:1_000
        @test sum(Pl(x,lmax=lmax)) >= 0
    end
end