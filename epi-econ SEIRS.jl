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
    S₁,E₁,I₁,R₁,S₂,E₂,I₂,R₂,Sᵤ,Eᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    α,β₁,β₂,βᵤ,γ,δ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂ = p

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ*S₁
    du[2] = dE₁ = β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Eᵤ - (δ+λ)*E₁
    du[3] = dI₁ = δ*E₁ - (γ+λ)*I₁
    du[4] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ*R₁
    du[5] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ*S₂
    du[6] = dE₂ = β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Eᵤ - (δ+λ)*E₂
    du[7] = dI₂ = δ*E₂ - (γ+λ)*I₂
    du[8] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ*R₂
    du[9] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ*(S₁+S₂)
    du[10] = dEᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - δ*Eᵤ - (P(Θ₁)+P(Θ₂))*Eᵤ + λ*(E₁+E₂)
    du[11] = dIᵤ = δ*Eᵤ - γ*Iᵤ + λ*(I₁+I₂)
    du[12] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ*(R₁+R₂)
    du[13] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(β₁-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ)/Q(Θ₁))/(c₁*η)
    du[14] = dΘ₂ = P(Θ₂)*((1-α)*(y₂-b-d*(β₂-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ)/Q(Θ₂))/(c₂*η)
end

tspan = (0,1000)
p = [0.5,0.4,0.3,0.2,1/7,1/3,0.0033,0.072,0.01,0.71,0.1,0.1,1,0.02,1.1,1.0]
#     α,  β₁, β₂, βᵤ, γ,  δ,    λ,    μ,   ρ,    b, c₁,c₂,d,r,y₁,y₂

# No infection
u₀ = [0.905, 0.0, 0.0, 0.0, 0.07, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 2.8, 0.0165]
prob = ODEProblem(myode, u₀, tspan, p)
sol = solve(prob,saveat=1)
Em0 = sol[1,:]+sol[2,:]+sol[3,:]
G0 = sol[4,:]+sol[5,:]+sol[6,:]
U0 = sol[7,:]+sol[8,:]+sol[9,:]

# Infection
u₀ = [0.905*0.99, 0.0, 0.905*0.01, 0.0, 0.07, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 2.8, 0.0165]
prob = ODEProblem(myode, u₀, tspan, p)
sol = solve(prob,saveat=1)
Em = sol[1,:]+sol[2,:]+sol[3,:]+sol[4,:]
G = sol[5,:]+sol[6,:]+sol[7,:]+sol[8,:]
U = sol[9,:]+sol[10,:]+sol[11,:]+sol[12,:]

S = sol[1,:]+sol[5,:]+sol[9,:]
E = sol[2,:]+sol[6,:]+sol[10,:]
I = sol[3,:]+sol[7,:]+sol[11,:]
R = sol[4,:]+sol[8,:]+sol[12,:]

# plot([E,G,U],ylims=(0,1))

# pe=plot([sol[1,:],sol[2,:],sol[3,:]])
# pg=plot([sol[4,:],sol[5,:],sol[6,:]])
# pu=plot([sol[7,:],sol[8,:],sol[9,:]])

# plot(pe,pg,pu,layout=(3,1))

pe=plot([Em,Em0],label=["Em" "Em0"])
pg=plot([G,G0],label=["G" "G0"])
pu=plot([U,U0],label=["U" "U0"])
psir=plot([S,E,I,R],label=["S" "E" "I" "R"])

plot(pe,pg,pu,psir,layout=(4,1))

# pout = plot(pe,pg,pu,psir,layout=(4,1))
# savefig(pout,"timeseries.pdf")
