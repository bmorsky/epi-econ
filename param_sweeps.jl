using DifferentialEquations, Plots, LaTeXStrings

# Output
E = zeros(101)
G = zeros(101)
U = zeros(101)
Θe = zeros(101)
Θg = zeros(101)

# Parameters
T = 2000 # final time

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*Θ^(1-η)
Q(Θ) = μ*Θ^(-η)

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ₁*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ #- (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ₁*R₁ + λ₂*R₂
    du[10] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(β₁-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁))/(c₁*η)
    du[11] = dΘ₂ = P(Θ₂)*((1-α)*(y₂+bg-b-d*(β₂-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂))/(c₂*η)
end
tspan = (0,T)

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]

# b sweep
for n=1:101
    b = 0.75*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b, c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pb = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.01, b, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]

    prob = ODEProblem(epiecon_ode, u₀, tspan, pb)
    sol = solve(prob,saveat=1)
    E[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θe[n] = sol[10,end]
    Θg[n] = sol[11,end]
end
sweep_b=plot(0:0.0075:0.75,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"b",legend=:outerright)

# β₂ sweep
for n=1:101
    β₂ = 0.25+0.1*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pbetag = [0.5, 0.4, β₂, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]

    prob = ODEProblem(epiecon_ode, u₀, tspan, pbetag)
    sol = solve(prob,saveat=1)
    E[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θe[n] = sol[10,end]
    Θg[n] = sol[11,end]
end
sweep_betag=plot(0.25:0.001:0.35,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_g",legend=:outerright)

# d sweep
for n=1:101
    d = (n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pd = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, d, 0.02, 1.1, 1.0, 0.0]

    prob = ODEProblem(epiecon_ode, u₀, tspan, pd)
    sol = solve(prob,saveat=1)
    E[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θe[n] = sol[10,end]
    Θg[n] = sol[11,end]
end
sweep_d=plot(0:0.01:1,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"d",legend=:outerright)

plot(sweep_betag,sweep_b,sweep_d,layout=(2,2))
savefig("sweep.pdf")
