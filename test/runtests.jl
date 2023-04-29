using TracyWidomBeta
using Test

@testset "TracyWidomBeta.jl" begin
    # Test typical inputs
    F_f=TW(2)
    F_s=TW(2;method="spectral",step="bdf5")
    TWG=TracyWidom
    @test F_f(-2) ≈ cdf(TWG,-2;beta=2,num_points=300) atol=1e-6
    @test F_s(-2) ≈ cdf(TWG,-2;beta=2,num_points=300) atol=1e-11

    # Test edge cases
    @test F_f(-10) ≈ 0 atol=1e-6
    @test F_s(-10) ≈ 0 atol=1e-11
    @test F_f(13/sqrt(2)) ≈ 1 atol=1e-6
    @test F_s(13/sqrt(2)) ≈ 1 atol=1e-11
end
