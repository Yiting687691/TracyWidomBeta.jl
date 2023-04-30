using TracyWidomBeta
using Test

@testset "TracyWidomBeta.jl" begin
    # Test typical inputs around the peak (cdf and pdf)
    F_f=TW(2)
    f_f=TW(2;pdf=true)
    F_s=TW(2;method="spectral",step="bdf5")
    f_s=TW(2;method="spectral",step="bdf5",pdf=true)
    @test F_f(-2) ≈ F_s(-2) atol=1e-6
    @test f_f(-2) ≈ f_s(-2) atol=1e-5

    # Test edge cases
    @test F_f(-10) ≈ F_s(-10) atol=1e-12
    @test F_f(13/sqrt(2)) ≈ F_s(13/sqrt(2)) atol=1e-12
    @test f_f(-10) ≈ f_s(-10) atol=1e-11
    @test f_f(13/sqrt(2)) ≈ f_s(13/sqrt(2)) atol=1e-11
end
