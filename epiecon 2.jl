using DifferentialEquations, Plots, LaTeXStrings

# Parameters
T = 600 # final time
#    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
p = [0.5, 0.4, 0.35, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.05, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]
p_gig = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.01]
p_quar = [0.5, 0.4, 0.3, 0.2, 1/7, 1.5*0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*Θ^(1-η)
Q(Θ) = μ*Θ^(-η)
W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

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
    du[10] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(I₁/(S₁+I₁+R₁)-Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁))/(c₁*η)
    du[11] = dΘ₂ = P(Θ₂)*((1-α)*(y₂+bg-b-d*(I₂/(S₂+I₂+R₂)-Iᵤ/(Sᵤ+Iᵤ+Rᵤ))) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂))/(c₂*η)
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

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]
prob = ODEProblem(epiecon_ode, u₀, tspan, p)
sol = solve(prob,saveat=1)
E = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = sol[7,:]+sol[8,:]+sol[9,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = sol[3,:]+sol[6,:]+sol[9,:]

Θe = sol[10,:]
Θg = sol[11,:]

pe=plot([E_DF,E],label=[L"E^{DF}" L"E"],xlabel=L"t",legend=:outerright)
pg=plot([G_DF,G],label=[L"G^{DF}" L"G"],xlabel=L"t",legend=:outerright)
pu=plot([U_DF,U],label=[L"U^{DF}" L"U"],xlabel=L"t",legend=:outerright)
psir=plot([S,I,R],label=[L"S" L"I" L"R"],xlabel=L"t",legend=:outerright)

ptightE=plot([Θe_DF,Θe],label=[L"\Theta_e^{DF}" L"\Theta_e"],xlabel=L"t",legend=:outerright)
ptightG=plot([Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel=L"t",legend=:outerright)

plot(pe,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("timeseries.pdf")

# Plot the effect of a government benefit to gig workers
prob = ODEProblem(epiecon_ode, u₀, tspan, p_gig)
sol = solve(prob,saveat=1)
E_gig = sol[1,:]+sol[2,:]+sol[3,:]
G_gig = sol[4,:]+sol[5,:]+sol[6,:]
U_gig = sol[7,:]+sol[8,:]+sol[9,:]

S_gig = sol[1,:]+sol[4,:]+sol[7,:]-S
I_gig = sol[2,:]+sol[5,:]+sol[8,:]-I
R_gig = sol[3,:]+sol[6,:]+sol[9,:]-R

Θe_gig = sol[10,:]
Θg_gig = sol[11,:]

pe_gig=plot([E_DF,E,E_gig],label=[L"E^{DF}" L"E" L"E^*"],xlabel=L"t",legend=:outerright)
pg_gig=plot([G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^*"],xlabel=L"t",legend=:outerright)
pu_gig=plot([U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^*"],xlabel=L"t",legend=:outerright)
psir_gig=plot([S_gig,I_gig,R_gig],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel=L"t",legend=:outerright)

ptightE_gig=plot([Θe_DF,Θe,Θe_gig],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
ptightG_gig=plot([Θg_DF,Θg,Θe_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)

plot(pe_gig,pg_gig,pu_gig,psir_gig,ptightE_gig,ptightG_gig,layout=(3,2))
savefig("timeseries_gig_benefits.pdf")

# Plot the effect of a government quarantine (λ₁>λ₂)
prob = ODEProblem(epiecon_ode, u₀, tspan, p_quar)
sol = solve(prob,saveat=1)
E_quar = sol[1,:]+sol[2,:]+sol[3,:]
G_quar = sol[4,:]+sol[5,:]+sol[6,:]
U_quar = sol[7,:]+sol[8,:]+sol[9,:]

S_quar = sol[1,:]+sol[4,:]+sol[7,:]-S
I_quar = sol[2,:]+sol[5,:]+sol[8,:]-I
R_quar = sol[3,:]+sol[6,:]+sol[9,:]-R

Θe_quar = sol[10,:]
Θg_quar = sol[11,:]

pe_quar=plot([E_DF,E,E_quar],label=[L"E^{DF}" L"E" L"E^*"],xlabel=L"t",legend=:outerright)
pg_quar=plot([G_DF,G,G_quar],label=[L"G^{DF}" L"G" L"G^*"],xlabel=L"t",legend=:outerright)
pu_quar=plot([U_DF,U,U_quar],label=[L"U^{DF}" L"U" L"U^*"],xlabel=L"t",legend=:outerright)
psir_quar=plot([S_quar,I_quar,R_quar],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel=L"t",legend=:outerright)

ptightE_quar=plot([Θe_DF,Θe,Θe_quar],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
ptightG_quar=plot([Θg_DF,Θg,Θe_quar],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)

plot(pe_quar,pg_quar,pu_quar,psir_quar,ptightE_quar,ptightG_quar,layout=(3,2))
savefig("timeseries_quarantine.pdf")

# function wage(yⱼ,βⱼ,I,Θ₁,Θ₂)
#     α,βᵤ,b,d,c₁,c₂ = 0.5,0.2,0.71,1,0.1,0.1
#     return α*yⱼ .+ (1-α)*(b .+ d*(βⱼ-βᵤ)*I) .+ α*c₁*Θ₁ .+ α*c₂*Θ₂
# end

# plot([wage(1.1,0.4,I,Θe,Θg),wage(1.0,0.3,I,Θe,Θg)],label=["We" "Wg"])
