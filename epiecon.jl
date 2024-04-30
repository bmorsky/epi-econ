using DifferentialEquations, Plots, LaTeXStrings

# Parameters
T = 600 # final time
p = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0]
#    α,   β₁,  β₂,  βᵤ,  γ,   λ,      μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*Θ^(1-η)
Q(Θ) = μ*Θ^(-η)

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂ = p

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ*(S₁+S₂)
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ*(I₁+I₂) #- (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ*(R₁+R₂)
    du[10] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(β₁-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ)/Q(Θ₁))/(c₁*η)
    du[11] = dΘ₂ = P(Θ₂)*((1-α)*(y₂-b-d*(β₂-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ)/Q(Θ₂))/(c₂*η)
end
tspan = (0,T)

# Disease free dynamics
u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0, 2.8, 0.0165]
prob = ODEProblem(epiecon_ode, u₀, tspan, p)
sol = solve(prob,saveat=1)
E_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = sol[7,:]+sol[8,:]+sol[9,:]
Θe_DF = sol[10,:]
Θg_DF = sol[11,:]

# 1% initially infectedInfection
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]
prob = ODEProblem(myode, u₀, tspan, p)
sol = solve(prob,saveat=1)
E = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = sol[7,:]+sol[8,:]+sol[9,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = sol[3,:]+sol[6,:]+sol[9,:]

Θe = sol[10,:]
Θg = sol[11,:]

# plot([E,G,U],ylims=(0,1))

# pe=plot([sol[1,:],sol[2,:],sol[3,:]])
# pg=plot([sol[4,:],sol[5,:],sol[6,:]])
# pu=plot([sol[7,:],sol[8,:],sol[9,:]])

# plot(pe,pg,pu,layout=(3,1))

pe=plot([E_DF,E],label=[L"E^{DF}" L"E"],xlabel=L"t",legend=:outerright)
pg=plot([G_DF,G],label=[L"G^{DF}" L"G"],xlabel=L"t",legend=:outerright)
pu=plot([U_DF,U],label=[L"U^{DF}" L"U"],xlabel=L"t",legend=:outerright)
psir=plot([S,I,R],label=[L"S" L"I" L"R"],xlabel=L"t",legend=:outerright)

ptightE=plot([Θe_DF,Θe],label=[L"\Theta_e^{DF}" L"\Theta_e"],xlabel=L"t",legend=:outerright)
ptightG=plot([Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel=L"t",legend=:outerright)

plot(pe,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("timeseries.pdf")

# plot(ptightE,ptightG,layout=(2,1))
# savefig("tightness.pdf")

# function wage(yⱼ,βⱼ,I,Θ₁,Θ₂)
#     α,βᵤ,b,d,c₁,c₂ = 0.5,0.2,0.71,1,0.1,0.1
#     return α*yⱼ .+ (1-α)*(b .+ d*(βⱼ-βᵤ)*I) .+ α*c₁*Θ₁ .+ α*c₂*Θ₂
# end

# plot([wage(1.1,0.4,I,Θe,Θg),wage(1.0,0.3,I,Θe,Θg)],label=["We" "Wg"])
