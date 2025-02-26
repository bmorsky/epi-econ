using DifferentialEquations, Groebner, LinearAlgebra, Symbolics, ModelingToolkit

#    α,   β₁,  β₂,  βᵤ,  γ,    λ,      μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
p = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]
α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p

# @independent_variables t
@variables  Vfes, Vfgs, Vfei, Vfgi, Vfer, Vfgr, PΘ₁, PΘ₂, QΘ₁, QΘ₂, I
# b, c₁, c₂, d, r, y₁, y₂, α, β₁, β₂, βᵤ, γ, λ, ρ,

Γs = PΘ₁*Vfes + PΘ₂*Vfgs
Γr = PΘ₁*Vfer + PΘ₂*Vfgr
VusVui = ((1-α)*d*(r+ρ) + α*((r+λ+ρ)*Γs - γ*Γr))/((1-α)*((r+γ)*(r+ρ) + (r+λ+ρ)*βᵤ*I))

Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I*(VusVui) + α*Γs
Wer = α*y₁ + (1 - α)*b + α*(PΘ₁*Vfer + PΘ₂*Vfgr)
Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I*(VusVui) + α*Γs
Wgr = α*y₂ + (1 - α)*b + α*(PΘ₁*Vfer + PΘ₂*Vfgr)
Wi = (1 − α)*b


eqns =  [y₁ - Wes - β₁*I*(Vfes - Vfei) - (r+λ)*Vfes,
         y₂ - Wgs - β₂*I*(Vfgs - Vfgi) - (r+λ)*Vfgs,
        -Wi - γ*(Vfei - Vfer) - (λ+r)*Vfei,
        -Wi - γ*(Vfgi - Vfgr) - (λ+r)*Vfgi,
        y₁ - Wer - ρ*(Vfer - Vfes) - (r+λ)*Vfer,
        y₂ - Wgr - ρ*(Vfgr - Vfgs) - (r+λ)*Vfgr]

sols = symbolic_linear_solve(eqns,[Vfes, Vfgs, Vfei, Vfgi, Vfer, Vfgr])