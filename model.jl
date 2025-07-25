using OrdinaryDiffEq, Parameters

# Parameters
T = 4000 # final time
#    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
# p = (α = 0.5, β₁= 0.6, β₂ = 0.5, βᵤ = 0.4,
#      γ = 1/7, λ₁ = 0.00046, λ₂ = 0.00046, μ = 0.01,
#      ρ = 0.01, b = 0.71, c₁ = 0.01, c₂ = 0.009, d = 1.0,
#      r = 0.00005, y₁ =  1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
     γ = 1/7, λ₁ = 0.0033, λ₂ = 0.0033, μ = 0.072,
     ρ = 0.01, b = 0.71, c₁ = 0.1, c₂ = 0.1, d = 0.1,
     r = 0.00002, y₁ =  1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
p_gig = merge(p,(;β₂ = 0.3, bg = 0.01)) # Copy p but change bg
# p_quar = merge(p,(;β₁= 0.9*0.4, y₁ = 0.99, λ₁ = 1.1*0.0033)) # Copy p but with changes to λ
p_quar = merge(p,(;β₂ = 0.3, λ₁ = 1.5*0.0033)) # Copy p but with changes to λ
# p_quar = merge(p,(;β₂ = 0.3, λ₁ = 1.5*0.00046)) # Copy p but with changes to λ

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    # Θ₂ = Θ₁
    @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η = p
    P(Θ) = μ*abs(Θ)^(1-η)
    Q(Θ) = μ*abs(Θ)^(-η)
    #W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

    # Θ₁ = (((1-α)*(y₁-b-d*(2*I₁/(S₁+I₁+R₁) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))))/(c₁*(r+λ₁)))^(1/η)
    # Θ₂ = (((1-α)*(y₂+bg-b-d*(2*I₂/(S₂+I₂+R₂) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) )/(c₂*(r+λ₂)))^(1/η)

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ₁*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ + P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ + P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ - (P(Θ₁)+P(Θ₂))*Iᵤ
    # du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ₁*R₁ + λ₂*R₂
    du[9] = 1.0-sum(u[1:9])
    # du[9] = dΘ = (1-α)*(y₁-b-d*(2*I₁/(S₁+I₁+R₁) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) - ((1-α)*(y₂+bg-b-d*(2*I₂/(S₂+I₂+R₂) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))))
    # du[10] = dΘ₁ = (1-α)*(y₁-b-d*(2*I₁/(S₁+I₁+R₁) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁)
    # du[11] = dΘ₂ = (1-α)*(y₂+bg-b-d*(2*I₂/(S₂+I₂+R₂) - Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂)
    du[10] = dΘ₁ = (1-α)*(y₁*((S₁+R₁)/(S₁+I₁+R₁)) - b) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁)
    du[11] = dΘ₂ = (1-α)*(y₂*((S₂+R₂)/(S₂+I₂+R₂)) - b) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂)
end
tspan = (0,T)
# M = zeros(11,11); [M[j,j] = 1 for j in [1:8...,10,11]]
# M = zeros(9,9); [M[j,j] = 1 for j in [1:8...]]
M = zeros(11,11); [M[j,j] = 1 for j in [1:8...]]
# Disease free dynamics
epiecon = ODEFunction(epiecon_ode,mass_matrix=M)
# epiecon = epiecon_ode