using ModelingToolkit

# Define units, parameters, and variables for ModelingToolkit
t = ModelingToolkit.t_nounits
D = ModelingToolkit.D_nounits
@parameters α,β₁,β₂,βᵤ,γ,η,λ₁,λ₂,μ,ρ,b₂,c₁,c₂,r,y₁,y₂
@variables begin
     S₁(..), [bounds=(0,1)]
     I₁(..), [bounds=(0,1)]
     R₁(..), [bounds=(0,1)]
     S₂(..), [bounds=(0,1)]
     I₂(..), [bounds=(0,1)]
     R₂(..), [bounds=(0,1)]
     Sᵤ(..), [bounds=(0,1)]
     Iᵤ(..), [bounds=(0,1)]
     Rᵤ(..), [bounds=(0,1)]
     Θ₁(..), [bounds=(0,100)]
     Θ₂(..), [bounds=(0,100)]
     J(..),
     bᵤ(..), [input = true, bounds = (0, 1)]
end

# Job loss and finding rates
P(Θ) = μ*(max(Θ, 1e-12))^(1-η)
Q(Θ) = μ*(max(Θ, 1e-12))^(-η)

# Equations, differential and algebraic
eqs = [   D(S₁(t)) ~ ρ*R₁(t) - β₁*S₁(t)*(I₁(t)+I₂(t)+Iᵤ(t)) + P(Θ₁(t))*Sᵤ(t) - λ₁*S₁(t),
          D(I₁(t)) ~ β₁*S₁(t)*(I₁(t)+I₂(t)+Iᵤ(t)) - (γ+λ₁)*I₁(t), #+ P(Θ₁)*Iᵤ
          D(R₁(t)) ~ γ*I₁(t) - ρ*R₁(t) + P(Θ₁(t))*Rᵤ(t) - λ₁*R₁(t),
          D(S₂(t)) ~ ρ*R₂(t) - β₂*S₂(t)*(I₁(t)+I₂(t)+Iᵤ(t)) + P(Θ₂(t))*Sᵤ(t) - λ₂*S₂(t),
          D(I₂(t)) ~ β₂*S₂(t)*(I₁(t)+I₂(t)+Iᵤ(t)) - (γ+λ₂)*I₂(t), #+ P(Θ₂)*Iᵤ
          D(R₂(t)) ~ γ*I₂(t) - ρ*R₂(t) + P(Θ₂(t))*Rᵤ(t) - λ₂*R₂(t),
          D(Sᵤ(t)) ~ ρ*Rᵤ(t) - βᵤ*Sᵤ(t)*(I₁(t)+I₂(t)+Iᵤ(t)) - (P(Θ₁(t))+P(Θ₂(t)))*Sᵤ(t) + λ₁*S₁(t)+λ₂*S₂(t),
          D(Iᵤ(t)) ~ βᵤ*Sᵤ(t)*(I₁(t)+I₂(t)+Iᵤ(t)) - γ*Iᵤ(t) + λ₁*I₁(t) + λ₂*I₂(t), # - (P(Θ₁)+P(Θ₂))*Iᵤ
          D(Rᵤ(t)) ~ γ*Iᵤ(t) - ρ*Rᵤ(t) - (P(Θ₁(t))+P(Θ₂(t)))*Rᵤ(t) + λ₁*R₁(t) + λ₂*R₂(t),
          0 ~ (1-α)*(y₁*((S₁(t)+R₁(t))/(S₁(t)+I₁(t)+R₁(t))) - bᵤ(t)) - α*c₁*Θ₁(t) - α*c₂*Θ₂(t) - c₁*(r+λ₁)/Q(Θ₁(t)),
          0 ~ (1-α)*(y₂*((S₂(t)+R₂(t))/(S₂(t)+I₂(t)+R₂(t))) + b₂ - bᵤ(t)) - α*c₁*Θ₁(t) - α*c₂*Θ₂(t) - c₂*(r+λ₂)/Q(Θ₂(t)),
          D(J(t)) ~ (I₁(t)+I₂(t)+Iᵤ(t))^2+(bᵤ(t)-0.4)^2
                    #0 ~ 1.0-S₁(t)-I₁(t)-R₁(t)-S₂(t)-I₂(t)-R₂(t)-Sᵤ(t)-Iᵤ(t)-Rᵤ(t),
     ]

# start and end times, costs, and constraints
(ts, te) = (0.0, 600.0)
costs = [J(te)]
cons = [bᵤ(te) ~ 0.4]#, S₁(t)+I₁(t)+R₁(t)+S₂(t)+I₂(t)+R₂(t)+Sᵤ(t)+Iᵤ(t)+Rᵤ(t) ~ 1]

# ModelingToolkit systems
@named epieconsys = System(eqs, t; costs, constraints = cons)
epieconsys = mtkcompile(epieconsys, inputs = [bᵤ(t)])

# Find initial conditions
include("model.jl")
# Disease-free scenario
u₀ = [SF₀, 0.0, 0.0, SG₀, 0.0, 0.0, SU₀, 0.0, 0.0, Θ₁₀, Θ₂₀]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = 1.0 .- F_DF .- G_DF
Θf_DF = sol[10,:]
Θg_DF = sol[11,:]
Wf_DF = W₁(0,Θf_DF,Θg_DF,p)
Wg_DF = W₂(0,Θf_DF,Θg_DF,p)
SW_DF = Wf_DF.*F_DF + p.bᵤ.*U_DF + (Wg_DF .+ p.b₂).*G_DF
# initial conditions for 1% initially infected
u0map = [ S₁(t) => SF₀*0.99, I₁(t) => SF₀*0.01, R₁(t) => 0.0, 
          S₂(t) => SG₀*0.99, I₂(t) => SG₀*0.01, R₂(t) => 0.0,
          Sᵤ(t) => SU₀*0.99, Iᵤ(t) => SU₀*0.01, Rᵤ(t) => 0.0,
          Θ₁(t) => Θ₁₀, Θ₂(t) => Θ₂₀, bᵤ(t) => 0.4, J(t) => 0.0]
# parameter values
pmap = [α => 0.5, β₁ => 0.4, β₂ => 0.38, βᵤ => 0.36, γ => 1/5, η => 0.5, λ₁ => 0.001096, λ₂ => 0.005479,
     μ => 0.01517, ρ => 0.01, b₂ => 0.2, c₁ => 0.17305498516997503, c₂ => 0.07253458196414889,
     r => 0.0001337, y₁ => 1.014065708418891, y₂ => 0.914065708418891]

# solve optimization problem
using InfiniteOpt, Ipopt, DiffEqDevTools
jprob = JuMPDynamicOptProblem(epieconsys, [u0map; pmap], (ts, te); dt = 1)
jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIIA5()));

# plot results
using Plots
pI = plot(jsol.sol,vars=((t, I₁, I₂, Iᵤ) -> (t, I₁+I₂+Iᵤ),0,2,5,8)) # plot total infections
pF = plot(jsol.sol,vars=((t, S₁, I₁, R₁) -> (t, S₁+I₁+R₁),0,7,8,9)) # plot total employmed in formal economy
pG = plot(jsol.sol,vars=((t, S₂, I₂, R₂) -> (t, S₂+I₂+R₂),0,4,6,6)) # plot total employmed in gig economy
pbu = plot(plot(jsol.input_sol)) # plot control, bᵤs