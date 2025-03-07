## Model

using NLsolve, OrdinaryDiffEq
using Parameters

p = (α = 0.5, β₁ = 0.4, β₂ = 0.3, βᵤ = 0.2, γ =  1/7,
     λ =  0.0033, μ =  0.072, ρ = 0.01, b = 0.71,
     c₁ = 0.1, c₂ = 0.1, d = 1.0, r = 0.02,
     y₁ = 1.5,y₂ = 1.4, bg = 0.0, η = 0.5)
p_gig = merge(p,(;bg = 0.01)) # Copy p but change bg
p_quar = merge(p,(;λ = 1.5*0.0033)) # Copy p but change λ

# Matching parameters and functions


# System of DAEs
function epiecon_ode!(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂,Vfes,Vfgs,Vfei,Vfgi,Vfer,Vfgr= u
    @unpack α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,η = p # Unpack macro to assign appropriate values to parameters

    Iₜ = I₁ + I₂ + Iᵤ
    τ = Sᵤ/(Sᵤ+Rᵤ)
    P(Θ) = μ*abs(Θ)^(1-η)
    Q(Θ) = μ*abs(Θ)^(-η)

    # vars = Ves, Vei, Ver, Vus, Vui, Vur, Vgs, Vgi, Vgr, Θ₁, Θ₂

    Γs = P(Θ₁)*Vfes + P(Θ₂)*Vfgs
    Γr = P(Θ₁)*Vfer + P(Θ₂)*Vfgr
    VusVui = ((1-α)*d*(r+ρ) + α*((r+γ+ρ)*Γs - γ*Γr))/((1-α)*((r+γ)*(r+ρ) + (r+γ+ρ)*βᵤ*Iₜ))

    # Wages
    Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*Iₜ*VusVui + α*Γs
    Wer = α*y₁ + (1 - α)*b + α*Γr
    Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*Iₜ*VusVui + α*Γs
    Wgr = α*y₂ + (1 - α)*b + α*Γr
    Wi = (1 − α)*b

    ## Differential equations
    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ*S₁+λ*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ*I₁ + λ*I₂ #- (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = dRᵤ = 1.0-S₁-S₂-Sᵤ-I₁-I₂-Iᵤ-R₁-R₂-Rᵤ # Conservation of population
    # du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ*R₁ + λ*R₂
    # Algebraic Equations
    #du[10] = dΘ₁ =  (r + λ)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*Iₜ*(VusVui) + (1-τ)*ρ*(Vfer-Vfes) + τ*β₁*Iₜ*(Vfes-Vfei)
    du[10] = dΘ₁ = Θ₁-(μ*(τ*Vfes+(1-τ)*Vfer)/c₁)^(1/η)
    #du[11] = dΘ₂ = (r + λ)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*Iₜ*(VusVui) + (1-τ)*ρ*(Vfgr-Vfgs) + τ*β₂*Iₜ*(Vfgs-Vfgi)
    du[11] =  dΘ₂ = Θ₂-(μ*(τ*Vfgs+(1-τ)*Vfgr)/c₂)^(1/η)
    du[12] = dVfes = y₁ - Wes - β₁*Iₜ*(Vfes - Vfei) - (r+λ)*Vfes
    du[13] = dVfgs = y₂ - Wgs - β₂*Iₜ*(Vfgs - Vfgi) - (r+λ)*Vfgs
    du[14] = dVfei = -Wi - γ*(Vfei - Vfer) - (λ+r)*Vfei
    du[15] = dVfgi = -Wi - γ*(Vfgi - Vfgr) - (λ+r)*Vfgi
    du[16] = dVfer = y₁ - Wer - ρ*(Vfer - Vfes) - (r+λ)*Vfer
    du[17] = dVfgr = y₂ - Wgr - ρ*(Vfgr - Vfgs) - (r+λ)*Vfgr
end

# Algebraic Equations
function equations!(F, vars, p,u)
    @unpack α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,η = p
    τ = u[7]/(u[7]+u[9])
    I₀ = sum(u[[2,5,8]])
    Θ₁,Θ₂,Vfes,Vfgs,Vfei,Vfgi,Vfer,Vfgr = vars

    P(Θ) = μ*abs(Θ)^(1-η)
    Q(Θ) = μ*abs(Θ)^(-η)

    # Vfes, Vfgs, Vfei, Vfgi, Vfer, Vfgr, Θ₁, Θ₂ = vars

    Γs = P(Θ₁)*Vfes + P(Θ₂)*Vfgs
    Γr = P(Θ₁)*Vfer + P(Θ₂)*Vfgr
    VusVui = ((1-α)*d*(r+ρ) + α*((r+γ+ρ)*Γs - γ*Γr))/((1-α)*((r+γ)*(r+ρ) + (r+γ+ρ)*βᵤ*I₀))

    # Wages
    Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I₀*(VusVui) + α*Γs
    Wer = α*y₁ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
    Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I₀*(VusVui) + α*Γs
    Wgr = α*y₂ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
    Wi = (1 − α)*b

    #F[1]] =   (r + λ)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*Iₜ*(VusVui) + (1-τ)*ρ*(Vfer-Vfes) + τ*β₁*Iₜ*(Vfes-Vfei)
    F[1] =  c₁-Q(Θ₁)*(τ*Vfes+(1-τ)*Vfer)
    #F[2] = (r + λ)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*Iₜ*(VusVui) + (1-τ)*ρ*(Vfgr-Vfgs) + τ*β₂*Iₜ*(Vfgs-Vfgi)
    F[2] =  dΘ₂ = c₂ - Q(Θ₂)*(τ*Vfgs+(1-τ)*Vfgr)
    F[3] = y₁ - Wes - β₁*I₀*(Vfes - Vfei) - (r+λ)*Vfes
    F[4] = y₂ - Wgs - β₂*I₀*(Vfgs - Vfgi) - (r+λ)*Vfgs
    F[5] = -Wi - γ*(Vfei - Vfer) - (λ+r)*Vfei
    F[6] = -Wi - γ*(Vfgi - Vfgr) - (λ+r)*Vfgi
    F[7] = y₁ - Wer - ρ*(Vfer - Vfes) - (r+λ)*Vfer
    F[8] = y₂ - Wgr - ρ*(Vfgr - Vfgs) - (r+λ)*Vfgr  
end




# Mass Matrix 
# DifferentialEquations.jl will solve systems of the type 
# M(du/dt) = f(u(t)) as a DAE
M = zeros(17,17)
for j in 1:8
    M[j,j] = 1.0
end

DAE = ODEFunction(epiecon_ode!,mass_matrix=M)
