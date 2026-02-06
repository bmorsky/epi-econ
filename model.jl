using DifferentialEquations, Parameters

# Parameters
T = 300 # final time
tspan = (0,T)

#    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
      γ = 1/7, λ₁ = 0.001096, λ₂ = 0.005479, μ = 0.01517,
      ρ = 0.01, b = 0.4, c₁ = 0.1731, c₂ = 0.07253, d = 1.0,
      r = 0.000137, y₁ =  1.0147, y₂ = 0.9141, bg = 0.2, η = 0.5,
      de = 0.0, dg = 0.0)
# p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
#      γ = 1/7, λ₁ = 0.001, λ₂ = 0.001, μ = 0.072,
#      ρ = 0.01, b = 0.2, c₁ = 0.2761, c₂ = 0.2761, d = 1.0,
#      r = 0.00013, y₁ =  1.139, y₂ = 1.139*0.95, bg = 0.0, η = 0.5)
# p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
#      γ = 1/7, λ₁ = 0.0033, λ₂ = 0.0033, μ = 0.072,
#      ρ = 0.01, b = 0.71, c₁ = 0.2761, c₂ = 0.2761, d = 1.0,
#      r = 0.004, y₁ =  1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
# p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
#      γ = 1/7, λ₁ = 0.0033, λ₂ = 0.0033, μ = 0.072,
#      ρ = 0.01, b = 0.71, c₁ = 0.1, c₂ = 0.1, d = 1.0,
#      r = 0.02, y₁ = 1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
# p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
#      γ = 1/7, λ₁ = 0.001, λ₂ = 0.001, μ = 0.01061,
#      ρ = 0.01, b = 0.71/7, c₁ = 0.04098, c₂ = 0.04098, d = 1.0,
#      r = 0.00013, y₁ =  1.014/7, y₂ = y₁*0.9, bg = 0.0, η = 0.5)
p_gig = merge(p,(;β₂ = 0.3, bg = 0.01)) # Copy p but change bg
p_quar = merge(p,(;β₂ = 0.3, λ₁ = 1.5*0.0033)) # Copy p but with changes to λ

function init_theta(du, u, p, t)
    ψ1,ψ2 = u
    @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η,de,dg = p
    Θ₁ = ψ1^(1/(1-η))
    Θ₂ = ψ2^(1/(1-η))
    du[1] = dΘ₁ = ((1-α)*(y₁ - b) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)*ψ1^(η/(1-η))/μ)
    du[2] = dΘ₂ = ((1-α)*(y₂-b+bg) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)*ψ2^(η/(1-η))/μ)
end
prob = SteadyStateProblem(init_theta, [1,1], p)
Ψ = solve(prob,SSRootfind())

function econ_ode(du, u, p, t)
    S₁,S₂,Sᵤ,ψ1,ψ2 = u
    @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η = p
    Θ₁ = ψ1^(1/(1-η))
    Θ₂ = ψ2^(1/(1-η))
    #W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

    du[1] = dS₁ = μ*ψ1*Sᵤ - λ₁*S₁
    du[2] = dS₂ = μ*ψ2*Sᵤ - λ₂*S₂
    du[3] = 1.0-sum(u[1:2])
    du[4] = dΘ₁ = ((1-α)*(y₁ - b) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)*ψ1^(η/(1-η))/μ)
    du[5] = dΘ₂ = ((1-α)*(y₂ - b+bg) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)*ψ2^(η/(1-η))/μ)
end
Mecon = zeros(5,5); [Mecon[j,j] = 1 for j in [1:2...]]
    # Disease free dynamics
econ = ODEFunction(econ_ode,mass_matrix=Mecon)
prob = ODEProblem(econ, [0.33,0.34,0.33,Ψ...], (0.0,2000.0), p)
sol = solve(prob,saveat=0.1, Rosenbrock23(), reltol=1e-8, abstol=1e-8)
SE₀ = sol[end][1]
SG₀ = sol[end][2]
SU₀ = 1.0 - SE₀ - SG₀
Θ1₀ = sol[end][4]
Θ2₀ = sol[end][5]

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,ψ1,ψ2 = u
    @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η,de,dg = p
    Θ₁ = ψ1^(1/(1-η))
    Θ₂ = ψ2^(1/(1-η))
    #W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) +μ*ψ1*Sᵤ - λ₁*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ +μ*ψ1*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + ψ1*Rᵤ*μ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + μ*ψ2*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ + μ*ψ2*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + μ*ψ2*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - μ*(ψ1+ψ2)*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ - μ*(ψ1+ψ2)*Iᵤ
    # du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ₁*R₁ + λ₂*R₂
    du[9] = 1.0-sum(u[1:9])
    du[10] = dΘ₁ = ((1-α)*(y₁*((S₁+R₁)/(S₁+I₁+R₁)) - b) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)*ψ1^(η/(1-η))/μ)
    du[11] = dΘ₂ = ((1-α)*(y₂*((S₂+R₂)/(S₂+I₂+R₂)) - b + bg) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)*ψ2^(η/(1-η))/μ)
end
M = zeros(11,11); [M[j,j] = 1 for j in [1:8...]]
# Disease free dynamics
epiecon = ODEFunction(epiecon_ode,mass_matrix=M)
