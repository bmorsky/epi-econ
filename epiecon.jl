## Make Figures 2,3,4

using Plots, LaTeXStrings
theme(:wong, lw = 2,
     size = (600,400),
     fontfamily = "Computer Modern",
     )

include("model.jl")

u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0, 2.8, 0.0165]
# u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1,abstol=1e-12)
E_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = 1.0 .- E_DF .- G_DF
Θe_DF = sol[10,:]
Θg_DF = sol[11,:]

# 1% initially infected
u₀ = [0.905*0.95, 0.905*0.05, 0.0, 0.07*0.95, 0.07*0.05, 0.0, 0.025*0.95, 0.025*0.05, 0.0, 2.8, 0.0165]
# u₀ = [0.905*0.95, 0.905*0.05, 0.0, 0.07*0.95, 0.07*0.05, 0.0, 0.025*0.95, 0.025*0.05, 0.0]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1)
E = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = 1.0 .- E .- G
# U = sol[7,:]+sol[8,:]+sol[9,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = 1.0 .- S .- I
# R = sol[3,:]+sol[6,:]+sol[9,:]

Θe = sol[10,:]
Θg = sol[11,:]
t = sol.t

pe=plot(t,[E_DF,E],label=[L"E^{DF}" L"E"],xlabel="",legend=:outerright)
pg=plot(t,[G_DF,G],label=[L"G^{DF}" L"G"],xlabel="",legend=:outerright)
pu=plot(t,[U_DF,U],label=[L"U^{DF}" L"U"],xlabel="",legend=:outerright)
psir=plot(t,[S,I,R],label=[L"S" L"I" L"R"],xlabel="",legend=:outerright)

ptightE=plot(t,[Θe_DF,Θe],label=[L"\Theta_e^{DF}" L"\Theta_e"],xlabel="time (days)",legend=:outerright)
ptightG=plot(t,[Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel="time (days)",legend=:outerright)

plot(pe,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("timeseries.pdf")

# Plot the effect of a government benefit to gig workers
prob = ODEProblem(epiecon, u₀, tspan, p_gig)
sol = solve(prob,saveat=0.1)
E_gig = sol[1,:]+sol[2,:]+sol[3,:]
G_gig = sol[4,:]+sol[5,:]+sol[6,:]
U_gig = 1.0 .- E_gig .- G_gig  #sol[7,:]+sol[8,:]+sol[9,:]

S_gig = sol[1,:]+sol[4,:]+sol[7,:]-S
I_gig = sol[2,:]+sol[5,:]+sol[8,:]-I
R_gig = 1.0 .- S_gig .- I_gig #sol[3,:]+sol[6,:]+sol[9,:]-R

Θe_gig = sol[10,:]
Θg_gig = sol[11,:]

pe_gig=plot(t,[E_DF,E,E_gig],label=[L"E^{DF}" L"E" L"E^*"],xlabel="",legend=:outerright)
pg_gig=plot(t,[G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^*"],xlabel="",legend=:outerright)
pu_gig=plot(t,[U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^*"],xlabel="",legend=:outerright)
psir_gig=plot(t,[S_gig,I_gig,R_gig],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel="",legend=:outerright)

ptightE_gig=plot(t,[Θe_DF,Θe,Θe_gig],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel="time (days)",legend=:outerright)
ptightG_gig=plot(t,[Θg_DF,Θg,Θg_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel="time (days)",legend=:outerright)



plot(pe_gig,pg_gig,pu_gig,psir_gig,ptightE_gig,ptightG_gig,layout=(3,2))
savefig("timeseries_gig_benefits.pdf")

# Plot the effect of a government quarantine (λ₁>λ₂)
prob = ODEProblem(epiecon, u₀, tspan, p_quar)
sol = solve(prob,saveat=0.1)
E_quar = sol[1,:]+sol[2,:]+sol[3,:]-E
G_quar = sol[4,:]+sol[5,:]+sol[6,:]-G
U_quar = 1.0 .- E_quar .- G_quar-U #sol[7,:]+sol[8,:]+sol[9,:]

S_quar = sol[1,:]+sol[4,:]+sol[7,:]-S
I_quar = sol[2,:]+sol[5,:]+sol[8,:]-I
R_quar = 1.0 .- S_quar .- I_quar -R

Θe_quar = sol[10,:]
Θg_quar = sol[11,:]

pe_quar=plot(t,[E_DF,E,E_quar],label=[L"E^{DF}" L"E" L"E^*"],xlabel="",legend=:outerright,)
pg_quar=plot(t,[G_DF,G,G_quar],label=[L"G^{DF}" L"G" L"G^*"],xlabel="",legend=:outerright)
pu_quar=plot(t,[U_DF,U,U_quar],label=[L"U^{DF}" L"U" L"U^*"],xlabel="",legend=:outerright)
psir_quar=plot(t,[S_quar,I_quar,R_quar],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel="",legend=:outerright)

ptightE_quar=plot(t,[Θe_DF,Θe,Θe_quar],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel="time (days)",legend=:outerright)
ptightG_quar=plot(t,[Θg_DF,Θg,Θg_quar],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel="time (days)",legend=:outerright)

plot(pe_quar,pg_quar,pu_quar,psir_quar,ptightE_quar,ptightG_quar,layout=(3,2))
savefig("timeseries_quarantine.pdf")

# function wage(yⱼ,βⱼ,I,Θ₁,Θ₂)
#     α,βᵤ,b,d,c₁,c₂ = 0.5,0.2,0.71,1,0.1,0.1
#     return α*yⱼ .+ (1-α)*(b .+ d*(βⱼ-βᵤ)*I) .+ α*c₁*Θ₁ .+ α*c₂*Θ₂
# end

# plot([wage(1.1,0.4,I,Θe,Θg),wage(1.0,0.3,I,Θe,Θg)],label=["We" "Wg"])

# p = (α = 0.5, β₁= 0.4, β₂ = 0.3, βᵤ = 0.2,
#      γ = 1/7, λ₁ = 0.0033, λ₂ = 0.0033, μ = 0.072,
#      ρ = 0.01, b = 0.71, c₁ = 0.1, c₂ = 0.1, d = 0.1,
#      r = 0.02, y₁ =  1.1, y₂ = 1.0, bg = 0.0, η = 0.5)
# function Θ₀(Θ,p)
#      @unpack α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg,η = p
#      return [(1-α)*(y₁-b) - α*c₁*Θ[1] - α*c₂*Θ[2] - c₁*(r+λ₁)/Q(Θ[1]),
#      (1-α)*(y₂+bg-b) - α*c₁*Θ[1] - α*c₂*Θ[2] - c₂*(r+λ₂)/Q(Θ[2])]
# end