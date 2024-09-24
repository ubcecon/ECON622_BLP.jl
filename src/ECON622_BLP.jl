module ECON622_BLP

import NLsolve

export share, delta, Integrate, share!

include("integrate.jl")

@doc raw"""
    share(δ, Σ, dFν, x)

Computes shares in random coefficient logit with mean tastes `δ`, observed characteristics `x`, unobserved taste distribution `dFν`, and taste covariances `Σ`. 
Assumes there is an outside option with u=0. The outside option has share `1-sum(s)`

# Arguments

- `δ` vector of length `J`
- `Σ` `K` by `K` matrix
- `x` `J` by `K` array
- `∫` AbstractIntegrator for integrating over distribution of `ν`

# Returns

- vector of length `J` consisting of $s_1$, ..., $s_J$
"""
function share(δ, Σ, x, ∫::Integrate.AbstractIntegrator)
  J,K = size(x)
  (length(δ) == J) || error("length(δ)=$(length(δ)) != size(x,1)=$J")
  (K,K) == size(Σ) || error("size(x,2)=$K != size(Σ)=$(size(Σ))")
  function shareν(ν)
    s = δ + x*Σ*ν
    smax=max(0,maximum(s))    
    return(exp.(s .- smax)./(sum(exp,s .- smax) + exp(-smax)))
  end
  return(∫(shareν))
end

function delta(s, Σ, x, ∫)
  δ = similar(s)
  delta!(δ, s, Σ, x, ∫)
  return(δ)
end
  
# to define custom reverse rule for enzyme
function delta!(δ, s, Σ, x, ∫)
  eqns! = let xi=x, Σi=Σ, ∫i=∫ 
    function(F,δ) 
        F .= s - share(δ,Σi,xi,∫i)
    end
  end
  δ0 = log.(s) .- log(1-sum(s))
  sol=NLsolve.nlsolve(eqns!, δ0, autodiff=:forward, method=:trust_region)
  if (sol.residual_norm > 1e-4)
    @warn "Possible problem in delta!(...)\n".*"$sol"
  end
  δ .= sol.zero
  return nothing
end 
  

import EnzymeCore 
import EnzymeCore.EnzymeRules
import EnzymeCore: Const, Duplicated, BatchDuplicated, Active, BatchDuplicatedNoNeed
import EnzymeCore.EnzymeRules: ConfigWidth, AugmentedReturn, needs_primal, needs_shadow
import Enzyme # needed for Enzyme.jacobian
import LinearAlgebra: lu

function EnzymeRules.forward(f::Const{typeof(delta)}, ::Type{<:Duplicated}, 
  s::Duplicated, Σ::Const, x::Const, ∫::Const)
  δ  = f.val(s.val, Σ.val, x.val, ∫.val)
  sf(d) = share(d, Σ.val, x.val, ∫.val)
  J = Enzyme.jacobian(Enzyme.Forward,sf, δ)
  dδ = J \ s.dval    
  return Duplicated(δ, dδ)
end

function EnzymeRules.forward(f::Const{typeof(delta)}, RT::Type{<:BatchDuplicated}, 
                             s::BatchDuplicated, Σ::Const, x::Const, ∫::Const)
  δ  = f.val(s.val, Σ.val, x.val, ∫.val)
  sf(d) = share(d, Σ.val, x.val, ∫.val)
  J = lu(Enzyme.jacobian(Enzyme.Forward,sf, δ))
  out = BatchDuplicated(δ, map(ds->J \ ds, s.dval))
  return out
end

function EnzymeRules.forward(f::Const{typeof(delta)}, RT::Type{<:BatchDuplicatedNoNeed}, 
                             s::BatchDuplicated{T,N}, Σ::Const, x::Const, ∫::Const) where {T, N}
  δ  = f.val(s.val, Σ.val, x.val, ∫.val)
  sf(d) = share(d, Σ.val, x.val, ∫.val)
  J = lu(Enzyme.jacobian(Enzyme.Forward,sf, δ))
  out = NTuple{N,T}(J \ ds for ds in s.dval)
  return out
end


function EnzymeRules.augmented_primal(config::ConfigWidth{W}, f::Const{typeof(delta!)}, RT::Type{<:Const}, 
                                      δ::Duplicated, s::Duplicated, Σ::Const, x::Const, ∫) where W
  f.val(δ.val,s.val, Σ.val, x.val, ∫.val)
  tape = (δ = copy(δ.val),s=copy(s.val),Σ=copy(Σ.val))
  @assert !needs_primal(config) && !needs_shadow(config)
  return AugmentedReturn(nothing , nothing, tape)
end

function EnzymeRules.reverse(config::ConfigWidth{W}, f::Const{typeof(delta!)}, dret::Type{<:Const}, tape, 
                             δ::Duplicated, s::Duplicated, Σ::Const, x::Const, ∫) where W
  sf(d) =share(d, tape.Σ, x.val, ∫.val)
  J = lu(Enzyme.jacobian(Enzyme.Forward,sf, tape.δ))
  s.dval .+= J \ δ.dval 
  return(nothing, nothing, nothing,nothing,nothing)
end

function broken!(x)
  out = 0
  for i ∈ eachindex(x)
    out+=exp(x[i]) 
  end
  x[1] = x[length(x)+1]
  return(out)
end



end