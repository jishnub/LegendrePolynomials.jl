using OffsetArrays,Test
using Legendre

@testset "Pl explicit" begin
	x = 2rand() - 1 
    P = Pl(x,lmax=5)
    @test P[0]≈1
    @test P[1]≈x
    @test P[2]≈(3x^2-1)/2
    @test P[3]≈(5x^3-3x)/2
    @test P[4]≈(35x^4-30x^2+3)/8
    @test P[5]≈(63x^5-70x^3+15x)/8
end

@testset "dPl explicit" begin
	x = 2rand() - 1 
    P = dPl(x,lmax=5)
    @test P[0]≈0
    @test P[1]≈1
    @test P[2]≈3x
    @test P[3]≈(-3+15x^2)/2
    @test P[4]≈1/8*(-60x + 140x^3)
    @test P[5]≈1/8*(15 - 210x^2 + 315x^4)
end

@testset "d2Pl explicit" begin
	x = 2rand() - 1 
    P = d2Pl(x,lmax=5)
    @test P[0]≈0
    @test P[1]≈0
    @test P[2]≈3
    @test P[3]≈15x
    @test P[4]≈1/8*(-60 + 420x^2)
    @test P[5]≈1/8*(-420x + 1260x^3)
end