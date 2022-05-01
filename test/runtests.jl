using Aqua
using Test
using LegendrePolynomials

@testset "Project quality assurance" begin
    Aqua.test_all(LegendrePolynomials, ambiguities = false)
end

using PerformanceTestTools

PerformanceTestTools.@include_foreach(
    "test_threaded.jl",
    [nothing,
    ["JULIA_NUM_THREADS" => Threads.nthreads() > 1 ? "1" : "2"],
    ],
)

