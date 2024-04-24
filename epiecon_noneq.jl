using DifferentialEquations, Plots

η = 0.5
μ = 0.072

function P(Θ)
    return μ*Θ^(1-η)
end

function Q(Θ)
    return μ*Θ^(-η)
end

function Qprime(Θ)
    return -η*μ*Θ^(-1-η)
end

function myode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂ = p

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ)*I₁
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ)*I₂
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ*(S₁+S₂)
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ*(I₁+I₂)
    du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ*(R₁+R₂)
    du[10] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(β₁-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ)/Q(Θ₁))/(c₁*η)
    du[11] = dΘ₂ = P(Θ₂)*((1-α)*(y₂-b-d*(β₂-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ)/Q(Θ₂))/(c₂*η)
end

tspan = (0,1000)
p = [0.5,0.4,0.3,0.2,1/7,0.0033,0.072,0.01,0.71,0.1,0.1,1,0.02,1.1,1.0]
#     α,  β₁, β₂, βᵤ, γ,    λ,    μ,   ρ,    b, c₁,c₂,d,r,y₁,y₂

# No infection
u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0, 2.8, 0.0165]
prob = ODEProblem(myode, u₀, tspan, p)
sol = solve(prob,saveat=1)
E0 = sol[1,:]+sol[2,:]+sol[3,:]
G0 = sol[4,:]+sol[5,:]+sol[6,:]
U0 = sol[7,:]+sol[8,:]+sol[9,:]
Θe0 = sol[10,:]
Θg0 = sol[11,:]

# Infection
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0, 2.8, 0.0165]
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

pe=plot([E,E0],label=["E" "E0"])
pg=plot([G,G0],label=["G" "G0"])
pu=plot([U,U0],label=["U" "U0"])
psir=plot([S,I,R],label=["S" "I" "R"])

ptightE=plot([Θe,Θe0],label=["Θe" "Θe0"])
ptightG=plot([Θg,Θg0],label=["Θg" "Θg0"])

plot(pe,pg,pu,psir,layout=(4,1))
savefig("timeseries.pdf")

plot(ptightE,ptightG,layout=(2,1))
savefig("tightness.pdf")