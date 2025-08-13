using OrdinaryDiffEq, Parameters

# Parameters
T = 2000 # final time
tspan = (0,T)
p = (α = 0.5, β₁= 0.4, β₂ = 0.38, βᵤ = 0.36, γ = 1/5, η = 0.5, λ₁ = 0.001096, λ₂ = 0.005479,
     μ = 0.01517, ρ = 0.01, b₂ = 0.2, bᵤ = 0.4, c₁ = 0.17305498516997503, c₂ = 0.07253458196414889,
     r = 0.000137, y₁ = 1.014065708418891, y₂ = 0.914065708418891)

p_gig = merge(p,(;b₂ = 0.21)) # Copy p but change b₂
p_jloss = merge(p,(;λ₁ = 1.5*0.001096)) # Copy p but with changes to λ₁

# Disease-free ODEs
function econ_ode(du, u, p, t)
    S₁,S₂,Sᵤ,Θ₁,Θ₂ = u
    @unpack α,β₁,β₂,βᵤ,γ,η,λ₁,λ₂,μ,ρ,b₂,bᵤ,c₁,c₂,r,y₁,y₂ = p
    P(Θ) = μ*(max(Θ, 1e-12))^(1-η)
    Q(Θ) = μ*(max(Θ, 1e-12))^(-η)

    du[1] = dS₁ = P(Θ₁)*Sᵤ - λ₁*S₁
    du[2] = dS₂ = P(Θ₂)*Sᵤ - λ₂*S₂
    du[3] = 1.0-sum(u[1:3])
    du[4] = dΘ₁ = (1-α)*(y₁ - bᵤ) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁)
    du[5] = dΘ₂ = (1-α)*(y₂ - bᵤ + b₂) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂)
end
Mecon = zeros(5,5); [Mecon[j,j] = 1 for j in [1:2...]]

# Solve disease-free dynamics to find initial conditions
econ = ODEFunction(econ_ode,mass_matrix=Mecon)
prob = ODEProblem(econ, [0.33,0.34,0.33,1.0,0.35341731768728646], (0.0,10000.0), p)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
SF₀ = sol[end][1]
SG₀ = sol[end][2]
SU₀ = 1.0 - SF₀ - SG₀
Θ₁₀ = sol[end][4]
Θ₂₀ = sol[end][5]

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    @unpack α,β₁,β₂,βᵤ,γ,η,λ₁,λ₂,μ,ρ,b₂,bᵤ,c₁,c₂,r,y₁,y₂ = p
    P(Θ) = μ*(max(Θ, 1e-12))^(1-η)
    Q(Θ) = μ*(max(Θ, 1e-12))^(-η)

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ₁*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ # - (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = 1.0-sum(u[1:9])
    du[10] = dΘ₁ = (1-α)*(y₁*((S₁+R₁)/(S₁+I₁+R₁)) - bᵤ) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁)
    du[11] = dΘ₂ = (1-α)*(y₂*((S₂+R₂)/(S₂+I₂+R₂)) + b₂ - bᵤ) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂)
end
M = zeros(11,11); [M[j,j] = 1 for j in [1:8...]]
epiecon = ODEFunction(epiecon_ode,mass_matrix=M)