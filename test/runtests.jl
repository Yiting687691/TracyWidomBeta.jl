using TracyWidomBeta
using Test

@testset "TracyWidomBeta.jl" begin
    # Test typical inputs
    F_f=TW(2)
    F_s=TW(2;method="spectral",step="bdf5")
    @test F_f(-2) ≈ F_s(-2) atol=1e-6

    # Test edge cases
    @test F_f(-10) ≈ F_s(-10) atol=1e-6
    @test F_f(13/sqrt(2)) ≈ F_s(13/sqrt(2)) atol=1e-6
end
