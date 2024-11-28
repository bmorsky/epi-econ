using DifferentialEquations, Plots, LaTeXStrings, NLsolve, DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra

# Parameters
T = 600 # final time
#    α,   β₁,  β₂,  βᵤ,  γ,   λ,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
p = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]
p_gig = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.01]
p_quar = [0.5, 0.4, 0.3, 0.2, 1/7, 1.5*0.0033, 0.0033, 0.072, 0.01, 0.71, 0.1, 0.1, 1, 0.02, 1.1, 1.0, 0.0]

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*abs(Θ)^(1-η)
Q(Θ) = μ*abs(Θ)^(-η)

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ = u
    α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p

    I = I₁ + I₂ + Iᵤ
    τ = Sᵤ/(Sᵤ+Rᵤ)

    # vars = Ves, Vei, Ver, Vus, Vui, Vur, Vgs, Vgi, Vgr, Θ₁, Θ₂

    # Define the system of equations
    function equations(F, vars)

        Ves, Vgs, Vei, Vgi, Ver, Vgr, Vus, Vui, Vur, Θ₁, Θ₂ = vars

        # Wages
        Vfer = (Ver-Vur)*(1 − α)/α
        Vfgr = (Vgr-Vur)*(1 − α)/α
        Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I*(Vus - Vui) + (1-α)*(P(Θ₁)*(Ves-Vus) + P(Θ₂)*(Vgs-Vus))
        Wer = α*y₁ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I*(Vus - Vui) + (1-α)*(P(Θ₁)*(Ves-Vus) + P(Θ₂)*(Vgs-Vus))
        Wgr = α*y₂ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wi = (1 − α)*b

        F[1] = Wes - β₁*I*(Ves - Vei) - λ*(Ves - Vus) - r*Ves
        F[2] = Wgs - β₂*I*(Vgs - Vei) - λ*(Vgs - Vus) - r*Vgs
        F[3] = Wi - d - γ*(Vei - Ver) - λ*(Vei - Vui) - r*Vei
        F[4] = Wi - d - γ*(Vgi - Vgr) - λ*(Vgi - Vui) - r*Vgi
        F[5] = Wer - ρ*(Ver - Ves) - λ*(Ver - Vur) - r*Ver
        F[6] = Wgr - ρ*(Vgr - Vgs) - λ*(Vgr - Vur) - r*Vgr
        F[7] = b - βᵤ*I*(Vus - Vui) - P(Θ₁)*(Vus - Ves) - P(Θ₂)*(Vus - Vgs) - r*Vus
        F[8] = b - d - γ*(Vui - Vur) - r*Vui
        F[9] = b - ρ*(Vur - Vus) - P(Θ₁)*(Vur - Ver) - P(Θ₂)*(Vur - Vgr) - r*Vur
        F[10] = ((r + λ)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*I*(Vus − Vui) + (1-α)*(1-τ)*ρ*(Ver-Vur − Ves+Vus)/α) + (1-α)*τ*β₁*I*(Ves-Vus − Vei+Vui)/α
        F[11] = ((r + λ)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*I*(Vus − Vui) + (1-α)*(1-τ)*ρ*(Vgr-Vur − Vgs+Vus)/α) + (1-α)*τ*β₂*I*(Vgs-Vus − Vgi+Vui)/α
    end

    # Initial guess for the variables
    initial_guess = [200, 200, 200, 200, 200, 200, 200, 200, 200, 1.0, 0.1]

    # Solve the system
    result = nlsolve(equations, initial_guess)

    Θ₁ = abs(result.zero[10])
    Θ₂ = abs(result.zero[11])

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ*S₁
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ*S₁+λ*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ*I₁ + λ*I₂ #- (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ*R₁ + λ*R₂
end
tspan = (0,T)

# Define a SavingCallback to store Θ₁ and Θ₂
# saved_values = SavedValues(Float64, Tuple{Float64, Float64}) # List to store results

function save_theta!(u, t, integrator)
    α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p
    S₁, I₁, R₁, S₂, I₂, R₂, Sᵤ, Iᵤ, Rᵤ = u

    I = I₁ + I₂ + Iᵤ
    τ = Sᵤ/(Sᵤ+Rᵤ)

    # Define the system of equations
    function equations(F, vars)

        Ves, Vgs, Vei, Vgi, Ver, Vgr, Vus, Vui, Vur, Θ₁, Θ₂ = vars

        # Wages
        Vfer = (Ver-Vur)*(1 − α)/α
        Vfgr = (Vgr-Vur)*(1 − α)/α
        Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I*(Vus - Vui) + (1-α)*(P(Θ₁)*(Ves-Vus) + P(Θ₂)*(Vgs-Vus))
        Wer = α*y₁ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I*(Vus - Vui) + (1-α)*(P(Θ₁)*(Ves-Vus) + P(Θ₂)*(Vgs-Vus))
        Wgr = α*y₂ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wi = (1 − α)*b

        F[1] = Wes - β₁*I*(Ves - Vei) - λ*(Ves - Vus) - r*Ves
        F[2] = Wgs - β₂*I*(Vgs - Vei) - λ*(Vgs - Vus) - r*Vgs
        F[3] = Wi - d - γ*(Vei - Ver) - λ*(Vei - Vui) - r*Vei
        F[4] = Wi - d - γ*(Vgi - Vgr) - λ*(Vgi - Vui) - r*Vgi
        F[5] = Wer - ρ*(Ver - Ves) - λ*(Ver - Vur) - r*Ver
        F[6] = Wgr - ρ*(Vgr - Vgs) - λ*(Vgr - Vur) - r*Vgr
        F[7] = b - βᵤ*I*(Vus - Vui) - P(Θ₁)*(Vus - Ves) - P(Θ₂)*(Vus - Vgs) - r*Vus
        F[8] = b - d - γ*(Vui - Vur) - r*Vui
        F[9] = b - ρ*(Vur - Vus) - P(Θ₁)*(Vur - Ver) - P(Θ₂)*(Vur - Vgr) - r*Vur
        F[10] = ((r + λ)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*I*(Vus − Vui) + (1-α)*(1-τ)*ρ*(Ver-Vur − Ves+Vus)/α) + (1-α)*τ*β₁*I*(Ves-Vus − Vei+Vui)/α
        F[11] = ((r + λ)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*I*(Vus − Vui) + (1-α)*(1-τ)*ρ*(Vgr-Vur − Vgs+Vus)/α) + (1-α)*τ*β₂*I*(Vgs-Vus − Vgi+Vui)/α
    end

    initial_guess = [200, 200, 200, 200, 200, 200, 200, 200, 200, 1.0, 0.1]
    result = nlsolve(equations, initial_guess)
    Θ₁ = abs(result.zero[10])
    Θ₂ = abs(result.zero[11])

    # Save time and Θ values
    return [Θ₁ Θ₂]
end

# Disease free dynamics
u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0]
saved_values = SavedValues(Float64, Array{Float64, 2})
save_cb = SavingCallback(save_theta!, saved_values, saveat=1:1:T)
prob = ODEProblem(epiecon_ode, u₀, tspan, p)
sol = solve(prob,saveat=1:1:T, callback=save_cb)
Θe_DF, Θg_DF = zeros(T), zeros(T)
for i=1:T
    Θe_DF[i] = saved_values.saveval[i][1]
    Θg_DF[i] = saved_values.saveval[i][2]
end

E_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = sol[7,:]+sol[8,:]+sol[9,:]

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0]
saved_values = SavedValues(Float64, Array{Float64, 2})
save_cb = SavingCallback(save_theta!, saved_values, saveat=1:1:T)
prob = ODEProblem(epiecon_ode, u₀, tspan, p)
sol = solve(prob, saveat=1:1:T, callback=save_cb)
Θe, Θg = zeros(T), zeros(T)
for i=1:T
    Θe[i] = saved_values.saveval[i][1]
    Θg[i] = saved_values.saveval[i][2]
end

E = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = sol[7,:]+sol[8,:]+sol[9,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = sol[3,:]+sol[6,:]+sol[9,:]

pe=plot([E_DF,E],label=[L"E^{DF}" L"E"],xlabel=L"t",legend=:outerright)
pg=plot([G_DF,G],label=[L"G^{DF}" L"G"],xlabel=L"t",legend=:outerright)
pu=plot([U_DF,U],label=[L"U^{DF}" L"U"],xlabel=L"t",legend=:outerright)
psir=plot([S,I,R],label=[L"S" L"I" L"R"],xlabel=L"t",legend=:outerright)

ptightE=plot([Θe_DF,Θe],label=[L"\Theta_e^{DF}" L"\Theta_e"],xlabel=L"t",legend=:outerright)
ptightG=plot([Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel=L"t",legend=:outerright)

plot(pe,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("timeseries.pdf")

# # Plot the effect of a government benefit to gig workers
# prob = ODEProblem(epiecon_ode, u₀, tspan, p_gig)
# sol = solve(prob,saveat=1)
# E_gig = sol[1,:]+sol[2,:]+sol[3,:]
# G_gig = sol[4,:]+sol[5,:]+sol[6,:]
# U_gig = sol[7,:]+sol[8,:]+sol[9,:]

# S_gig = sol[1,:]+sol[4,:]+sol[7,:]-S
# I_gig = sol[2,:]+sol[5,:]+sol[8,:]-I
# R_gig = sol[3,:]+sol[6,:]+sol[9,:]-R

# Θe_gig = sol[10,:]
# Θg_gig = sol[11,:]

# pe_gig=plot([E_DF,E,E_gig],label=[L"E^{DF}" L"E" L"E^*"],xlabel=L"t",legend=:outerright)
# pg_gig=plot([G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^*"],xlabel=L"t",legend=:outerright)
# pu_gig=plot([U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^*"],xlabel=L"t",legend=:outerright)
# psir_gig=plot([S_gig,I_gig,R_gig],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel=L"t",legend=:outerright)

# ptightE_gig=plot([Θe_DF,Θe,Θe_gig],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
# ptightG_gig=plot([Θg_DF,Θg,Θe_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)

# plot(pe_gig,pg_gig,pu_gig,psir_gig,ptightE_gig,ptightG_gig,layout=(3,2))
# savefig("timeseries_gig_benefits.pdf")

# # Plot the effect of a government quarantine (λ₁>λ₂)
# prob = ODEProblem(epiecon_ode, u₀, tspan, p_quar)
# sol = solve(prob,saveat=1)
# E_quar = sol[1,:]+sol[2,:]+sol[3,:]
# G_quar = sol[4,:]+sol[5,:]+sol[6,:]
# U_quar = sol[7,:]+sol[8,:]+sol[9,:]

# S_quar = sol[1,:]+sol[4,:]+sol[7,:]-S
# I_quar = sol[2,:]+sol[5,:]+sol[8,:]-I
# R_quar = sol[3,:]+sol[6,:]+sol[9,:]-R

# Θe_quar = sol[10,:]
# Θg_quar = sol[11,:]

# pe_quar=plot([E_DF,E,E_quar],label=[L"E^{DF}" L"E" L"E^*"],xlabel=L"t",legend=:outerright)
# pg_quar=plot([G_DF,G,G_quar],label=[L"G^{DF}" L"G" L"G^*"],xlabel=L"t",legend=:outerright)
# pu_quar=plot([U_DF,U,U_quar],label=[L"U^{DF}" L"U" L"U^*"],xlabel=L"t",legend=:outerright)
# psir_quar=plot([S_quar,I_quar,R_quar],label=[L"\Delta S" L"\Delta I" L"\Delta R"],xlabel=L"t",legend=:outerright)

# ptightE_quar=plot([Θe_DF,Θe,Θe_quar],label=[L"\Theta_e^{DF}" L"\Theta_e" L"\Theta_e^*"],xlabel=L"t",legend=:outerright)
# ptightG_quar=plot([Θg_DF,Θg,Θe_quar],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^*"],xlabel=L"t",legend=:outerright)

# plot(pe_quar,pg_quar,pu_quar,psir_quar,ptightE_quar,ptightG_quar,layout=(3,2))
# savefig("timeseries_quarantine.pdf")

# function wage(yⱼ,βⱼ,I,Θ₁,Θ₂)
#     α,βᵤ,b,d,c₁,c₂ = 0.5,0.2,0.71,1,0.1,0.1
#     return α*yⱼ .+ (1-α)*(b .+ d*(βⱼ-βᵤ)*I) .+ α*c₁*Θ₁ .+ α*c₂*Θ₂
# end

# plot([wage(1.1,0.4,I,Θe,Θg),wage(1.0,0.3,I,Θe,Θg)],label=["We" "Wg"])
