module Integrate

using Distributions
import Sobol: skip, SobolSeq
import Base.Iterators: take, Repeated, product, repeated
import HCubature: hcubature
import LinearAlgebra: cholesky
using FastGaussQuadrature, LinearAlgebra
using SparseGrids


abstract type AbstractIntegrator end

(∫::AbstractIntegrator)(f::Function) = sum(w*f(x) for (w,x) in zip(∫.w, ∫.x))

struct FixedNodeIntegrator{Tx,Tw} <: AbstractIntegrator
    x::Tx
    w::Tw
end

MonteCarloIntegrator(distribution::Distribution, ndraw=100)=FixedNodeIntegrator([rand(distribution) for i=1:ndraw], Repeated(1/ndraw))

function QuasiMonteCarloIntegrator(distribution::UnivariateDistribution, ndraws=100)
    ss = skip(SobolSeq(1), ndraw)
    x = [quantile(distribution, x) for x in take(ss,ndraw)]
    w = Repeated(1/ndraw)
    FixedNodeIntegrator(x,w)
end 

function QuasiMonteCarloIntegrator(distribution::AbstractMvNormal, ndraw=100)
    ss = skip(SobolSeq(length(distribution)), ndraw)
    L = cholesky(distribution.Σ).L
    x = [L*quantile.(Normal(), x) for x in take(ss,ndraw)]
    w = Repeated(1/ndraw)
    FixedNodeIntegrator(x,w)
end 

function  QuadratureIntegrator(dx::AbstractMvNormal; ndraw=100)
  n = Int(ceil(ndraw^(1/length(dx))))
  x, w = gausshermite(n)
  L = cholesky(dx.Σ).L
  x = vec([√2*L*vcat(xs...) + dx.μ for xs in product(repeated(x,length(dx))...)])
  w = vec([prod(ws)/π^(length(dx)/2) for ws in product(repeated(w, length(dx))...)])
  FixedNodeIntegrator(x,w)
end

function SparseGridIntegrator(dx::AbstractMvNormal; order=5)
  X, W = sparsegrid(length(dx), order, gausshermite, sym=true)
  L = cholesky(dx.Σ).L
  X = [√2*L*x + dx.μ for x ∈ X]
  W /= π^(length(dx)/2)
  FixedNodeIntegrator(X,W)
end

struct AdaptiveIntegrator{FE,FT,FJ,A,L} <: AbstractIntegrator
    eval::FE
    transform::FT
    detJ::FJ
    args::A
    limits::L
end

(∫::AdaptiveIntegrator)(f::Function) = ∫.eval(t->f(∫.transform(t))*∫.detJ(t), ∫.limits...; ∫.args...)[1]

function AdaptiveIntegrator(dist::AbstractMvNormal; eval=hcubature, options=())
    D = length(dist)
    x(t) = t./(1 .- t.^2)
    Dx(t) = prod((1 .+ t.^2)./(1 .- t.^2).^2)*pdf(dist,x(t))
    args = options
    limits = (-ones(D), ones(D))
    AdaptiveIntegrator(hcubature,x,Dx,args, limits)
end

end