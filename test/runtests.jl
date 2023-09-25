using Test, TestItemRunner


#=
 # Just do this once
import Pkg
Pkg.activate()
Pkg.develop(path=normpath(joinpath(@__DIR__, "..")))
Pkg.activate(@__DIR__)
 =#


@testitem "integrator" begin
    import ECON622_BLP as BLP
    import ECON622_BLP: Integrate
    using Distributions, LinearAlgebra
    # includet("../src/integrate.jl") # for interactive execution
    
    dimx = 3
    A = rand(dimx,dimx)
    Σ = A*A' + I
    dx = MvNormal(zeros(dimx), Σ)
    ∫a = Integrate.AdaptiveIntegrator(dx, options=(rtol=1e-4,initdiv=3))        
    V = ∫a(x->x*x')
    @test isapprox(V, Σ, rtol=1e-5)
    
    f(x) = exp(x[1])/sum(exp.(x)) 
    val = ∫a(f)
    for N ∈ [1_000, 10_000]
        ∫mc = Integrate.MonteCarloIntegrator(dx, N)
        ∫qmc = Integrate.QuasiMonteCarloIntegrator(dx,N)        
        @test isapprox(∫mc(x->x*x'), Σ, rtol=10/sqrt(N))
        @test isapprox(∫qmc(x->x*x'), Σ, rtol=10/sqrt(N))
        
        @test isapprox(∫mc(f),val,rtol=2/sqrt(N))
        @test isapprox(∫qmc(f),val,rtol=2/sqrt(N))
    end
    for N ∈ [25, 100, 200]
        ∫ = Integrate.QuadratureIntegrator(dx,ndraw=N)
        n=ceil(N^(1/dimx))
        @test isapprox(∫(x->x*x'), Σ, rtol=1/2^n)
        @test isapprox(∫(f),val,rtol=1/2^n)
    end

    for order ∈ [3, 5, 8, 10]
        ∫ = Integrate.SparseGridIntegrator(dx,order=order)
        n=length(∫.w)
        @test isapprox(∫(x->x*x'), Σ, rtol=1/n)
        @test isapprox(∫(f),val,rtol=log(n)/n)
    end
end

@testitem "share=δ⁻¹" begin    
    import ECON622_BLP as BLP
    import ECON622_BLP: Integrate
    using Distributions, LinearAlgebra

    J = 4
    dimx=2
    dx = MvNormal(dimx, 1.0)
    Σ = [1 0.5; 0.5 1]
    N = 50
    ∫ = Integrate.QuasiMonteCarloIntegrator(dx, N)
    X = [(-1.).^(1:J) 1:J]
    δ = collect((1:J)./J)
    s = BLP.share(δ,Σ,X,∫) 
    d = BLP.delta(s, Σ, X, ∫)
    @test d ≈ δ

    J = 10
    dimx = 4
    X = rand(J, dimx)
    dx = MvNormal(dimx, 1.0)
    Σ = I + ones(dimx,dimx)
    ∫ = Integrate.QuasiMonteCarloIntegrator(dx, N)
    δ = 1*rand(J)
    s = BLP.share(δ,Σ,X,∫) 
    d = BLP.delta(s, Σ, X, ∫)
    @test isapprox(d, δ, rtol=1e-6)    

    @testset "hard delta" begin 
        δ = 5*rand(J) 
        @show s = BLP.share(δ,Σ,X,∫)
        d = BLP.delta(s, Σ, X, ∫)
        @test_skip isapprox(d, δ, rtol=1e-4)    
    end
end

@testitem "differentialibility" begin
    import ECON622_BLP as BLP
    import ECON622_BLP: Integrate
    using Distributions, LinearAlgebra

    import FiniteDiff, ForwardDiff, Enzyme, Zygote
    import Enzyme: Active, Const, Duplicated

    J = 4
    dimx=2
    dx = MvNormal(dimx, 1.0)
    Σ = [1 0.5; 0.5 1]
    ∫ = Integrate.QuasiMonteCarloIntegrator(dx, 40)
    X = [(-1.).^(1:J) 1:J]
    δ = collect((1:J)./J)
    s = BLP.share(δ, Σ, X, ∫)
    sharefn = d->BLP.share(d,Σ,X,∫)
    @testset "share function" begin
        Jfd = FiniteDiff.finite_difference_jacobian(sharefn, δ)
        @test isapprox(Jfd, ForwardDiff.jacobian(sharefn, δ), rtol=1e-4)
        @test isapprox(Jfd, Zygote.jacobian(sharefn, δ)[1], rtol=1e-4)
        Enzyme.API.runtimeActivity!(true)
        Jer = Enzyme.jacobian(Enzyme.Reverse, sharefn, δ, Val(J))
        @test isapprox(Jfd, Jer, rtol=1e-4)
        Jef = Enzyme.jacobian(Enzyme.Forward, sharefn, δ, Val(J))
        @test isapprox(Jfd, Jef, rtol=1e-4)
    end 

    @testset "delta function" begin
        dfn = s->BLP.delta(s,Σ,X,∫)
        Jfd = FiniteDiff.finite_difference_jacobian(dfn, s)
        @test isapprox(Jfd, ForwardDiff.jacobian(dfn, s), rtol=1e-4)
        @test isapprox(Jfd, Enzyme.jacobian(Enzyme.Forward, dfn, s), rtol=1e-4)
        @test isapprox(Jfd, Enzyme.jacobian(Enzyme.Reverse, dfn, s, Val(J)), rtol=1e-4)
    end
  
end

@run_package_tests
