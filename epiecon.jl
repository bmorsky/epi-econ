using OrdinaryDiffEq, Plots, LaTeXStrings
using Parameters

theme(:bright, lw = 2,
     size = (600,320),
     fontfamily = "Computer Modern")


# Parameters
T = 600 # final time
#    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
     γ = 1/7, λ₁ = 0.0033, λ₂ = 0.0033, μ = 0.072,
     ρ = 0.01, b = 0.71, c₁ = 0.1, c₂ = 0.1, d = 1.0,
     r = 0.02, y₁ =  1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
p_gig = merge(p,(;β₂ = 0.3, bg = 0.01)) # Copy p but change bg
p_quar = merge(p,(;β₂ = 0.3, λ₁ = 1.5*0.0033)) # Copy p but with changes to λ


# Matching parameters and functions

#P(Θ) = μ*Θ^(1-η)
#Q(Θ) = μ*Θ^(-η)
#W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ,Θ₁,Θ₂ = u
    @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η = p
    P(Θ) = μ*Θ^(1-η)
    Q(Θ) = μ*Θ^(-η)
    #W(b,dI) = V_{j,s}\dot{S}_j + V_{j,i}\dot{I}_j + V_{j,r}\dot{R}_j - V_{u,s}\dot{S}_u - V_{u,i}\dot{I}_u - V_{u,r}\dot{R}_u

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ₁*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ #- (P(Θ₁)+P(Θ₂))*Iᵤ
    #du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ₁*R₁ + λ₂*R₂
    du[9] = 1.0-sum(u[1:9])
    du[10] = dΘ₁ = P(Θ₁)*((1-α)*(y₁-b-d*(β₁-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₁*(r+λ₁)/Q(Θ₁))/(c₁*η)
    du[11] = dΘ₂ = P(Θ₂)*((1-α)*(y₂+bg-b-d*(β₂-βᵤ)*(I₁+I₂+Iᵤ)) - α*c₁*Θ₁ - α*c₂*Θ₂ - c₂*(r+λ₂)/Q(Θ₂))/(c₂*η)
end
tspan = (0,T)
M = zeros(11,11); [M[j,j] = 1 for j in [1:8...,10,11]]
# Disease free dynamics
u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0, 2.8, 0.0165]
epiecon = ODEFunction(epiecon_ode,mass_matrix=M)
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1,abstol=1e-12)
E_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = sol[7,:]+sol[8,:]+sol[9,:]
Θe_DF = sol[10,:]
Θg_DF = sol[11,:]

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1,abstol=1e-12,reltol=1e-8)
E = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = sol[7,:]+sol[8,:]+sol[9,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = sol[3,:]+sol[6,:]+sol[9,:]

Θe = sol[10,:]
Θg = sol[11,:]
t = sol.t

pe=plot(t,[E_DF,E],label=[L"E^{DF}" L"E"],xlabel="",legend=:outerright)
pg=plot(t,[G_DF,G],label=[L"G^{DF}" L"G"],xlabel="",legend=:outerright)
pu=plot(t,[U_DF,U],label=[L"U^{DF}" L"U"],xlabel="",legend=:outerright)
psir=plot(t,[S,I,R],label=[L"S" L"I" L"R"],xlabel="",legend=:outerright)

ptightE=plot(t,[Θe_DF,Θe],label=[L"\Theta_e^{DF}" L"\Theta_e"],xlabel=L"t",legend=:outerright)
ptightG=plot(t,[Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel=L"t",legend=:outerright)

plot(pe,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("Figures/timeseries.pdf")

# Plot the effect of a government benefit to gig workers
prob = ODEProblem(epiecon, u₀, tspan, p_gig)
sol = solve(prob,saveat=0.1,reltol=1e-8)
E_gig = sol[1,:]+sol[2,:]+sol[3,:]
G_gig = sol[4,:]+sol[5,:]+sol[6,:]
U_gig = sol[7,:]+sol[8,:]+sol[9,:]

S_gig = sol[1,:]+sol[4,:]+sol[7,:]-S
I_gig = sol[2,:]+sol[5,:]+sol[8,:]-I
R_gig = sol[3,:]+sol[6,:]+sol[9,:]-R

Θe_gig = sol[10,:]
Θg_gig = sol[11,:]

pe_gig=plot(t,[E_DF,E,E_gig],label=[L"E^{DF}" L"E" L"E^*"],xlabel="",legend=:outerright)
pg_gig=plot(t,[G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^*"],xlabel="",legend=:outerright)
pu_gig=plot(t,[U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^*"],xlabel="",legend=:outerright)
psir_gig=plot(t,[S_gig,I_gig,R_gig],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel="",legend=:outerright)

ptightE_gig=plot(t,[Θe_DF,Θe,Θe_gig],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
ptightG_gig=plot(t,[Θg_DF,Θg,Θg_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)



plot(pe_gig,pg_gig,pu_gig,psir_gig,ptightE_gig,ptightG_gig,layout=(3,2))
savefig("Figures/timeseries_gig_benefits.pdf")

# Plot the effect of a government quarantine (λ₁>λ₂)
prob = ODEProblem(epiecon, u₀, tspan, p_quar)
sol = solve(prob,saveat=0.1,reltol=1e-8)
E_quar = sol[1,:]+sol[2,:]+sol[3,:]
G_quar = sol[4,:]+sol[5,:]+sol[6,:]
U_quar = sol[7,:]+sol[8,:]+sol[9,:]

S_quar = sol[1,:]+sol[4,:]+sol[7,:]-S
I_quar = sol[2,:]+sol[5,:]+sol[8,:]-I
R_quar = sol[3,:]+sol[6,:]+sol[9,:]-R

Θe_quar = sol[10,:]
Θg_quar = sol[11,:]

pe_quar=plot(t,[E_DF,E,E_quar],label=[L"E^{DF}" L"E" L"E^*"],xlabel="",legend=:outerright,)
pg_quar=plot(t,[G_DF,G,G_quar],label=[L"G^{DF}" L"G" L"G^*"],xlabel="",legend=:outerright)
pu_quar=plot(t,[U_DF,U,U_quar],label=[L"U^{DF}" L"U" L"U^*"],xlabel="",legend=:outerright)
psir_quar=plot(t,[S_quar,I_quar,R_quar],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel="",legend=:outerright)

ptightE_quar=plot(t,[Θe_DF,Θe,Θe_quar],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
ptightG_quar=plot(t,[Θg_DF,Θg,Θg_quar],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)

plot(pe_quar,pg_quar,pu_quar,psir_quar,ptightE_quar,ptightG_quar,layout=(3,2))
savefig("Figures/timeseries_quarantine.pdf")

# function wage(yⱼ,βⱼ,I,Θ₁,Θ₂)
#     α,βᵤ,b,d,c₁,c₂ = 0.5,0.2,0.71,1,0.1,0.1
#     return α*yⱼ .+ (1-α)*(b .+ d*(βⱼ-βᵤ)*I) .+ α*c₁*Θ₁ .+ α*c₂*Θ₂
# end

# plot([wage(1.1,0.4,I,Θe,Θg),wage(1.0,0.3,I,Θe,Θg)],label=["We" "Wg"])
