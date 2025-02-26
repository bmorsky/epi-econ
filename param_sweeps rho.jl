using DifferentialEquations, Plots, LaTeXStrings, NLsolve, DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra

# Output
E = zeros(11)
G = zeros(11)
U = zeros(11)
Θe = zeros(11)
Θg = zeros(11)

# Parameters
T = 2000 # final time

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*abs(Θ)^(1-η)
Q(Θ) = μ*abs(Θ)^(-η)

# System of ODEs
function epiecon_ode(du, u, p, t)
    S₁,I₁,R₁,S₂,I₂,R₂,Sᵤ,Iᵤ,Rᵤ = u
    α,β₁,β₂,βᵤ,γ,λ₁,λ₂,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p

    I = I₁ + I₂ + Iᵤ
    τ = Sᵤ/(Sᵤ+Rᵤ)

    # vars = Ves, Vei, Ver, Vus, Vui, Vur, Vgs, Vgi, Vgr, Θ₁, Θ₂

    # Define the system of equations
    function equations(F, vars)

        Ves, Vgs, Vei, Vgi, Ver, Vgr, Vus, Vui, Vur, Θ₁, Θ₂ = vars

        # Wages
        Vfes = (Ves-Vus)*(1 − α)/α
        Vfgs = (Vgs-Vus)*(1 − α)/α
        Vfei = (Vei-Vui)*(1 − α)/α
        Vfgi = (Vgi-Vui)*(1 − α)/α
        Vfer = (Ver-Vur)*(1 − α)/α
        Vfgr = (Vgr-Vur)*(1 − α)/α
        Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I*(Vus - Vui) + α*(P(Θ₁)*Vfes + P(Θ₂)*Vfgs)
        Wer = α*y₁ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I*(Vus - Vui) + α*(P(Θ₁)*Vfes + P(Θ₂)*Vfgs)
        Wgr = α*y₂ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
        Wi = (1 − α)*b

        F[1] = Wes - β₁*I*(Ves - Vei) - λ₁*(Ves - Vus) - r*Ves
        F[2] = Wgs - β₂*I*(Vgs - Vei) - λ₂*(Vgs - Vus) - r*Vgs
        F[3] = Wi - d - γ*(Vei - Ver) - λ₁*(Vei - Vui) - r*Vei
        F[4] = Wi - d - γ*(Vgi - Vgr) - λ₂*(Vgi - Vui) - r*Vgi
        F[5] = Wer - ρ*(Ver - Ves) - λ₁*(Ver - Vur) - r*Ver
        F[6] = Wgr - ρ*(Vgr - Vgs) - λ₂*(Vgr - Vur) - r*Vgr
        F[7] = b - βᵤ*I*(Vus - Vui) - P(Θ₁)*(Vus - Ves) - P(Θ₂)*(Vus - Vgs) - r*Vus
        F[8] = b - d - γ*(Vui - Vur) - r*Vui
        F[9] = b - ρ*(Vur - Vus) - P(Θ₁)*(Vur - Ver) - P(Θ₂)*(Vur - Vgr) - r*Vur
        F[10] = (r + λ₁)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*I*(Vus − Vui) + (1-τ)*ρ*(Vfer-Vfes) + τ*β₁*I*(Vfes-Vfei)
        F[11] = (r + λ₂)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*I*(Vus − Vui) + (1-τ)*ρ*(Vfgr-Vfgs) + τ*β₂*I*(Vfgs-Vfgi)
    end

    # Initial guess for the variables
    initial_guess = [300, 300, 300, 300, 300, 300, 300, 300, 300, 2.0, 0.01]

    # Solve the system
    result = nlsolve(equations, initial_guess)

    Θ₁ = result.zero[10]
    Θ₂ = result.zero[11]

    du[1] = dS₁ = ρ*R₁ - β₁*S₁*(I₁+I₂+Iᵤ) + P(Θ₁)*Sᵤ - λ₁*S₁ #(1+cos(2*pi*t/1000))
    du[2] = dI₁ = β₁*S₁*(I₁+I₂+Iᵤ) - (γ+λ₁)*I₁ #+ P(Θ₁)*Iᵤ
    du[3] = dR₁ = γ*I₁ - ρ*R₁ + P(Θ₁)*Rᵤ - λ₁*R₁
    du[4] = dS₂ = ρ*R₂ - β₂*S₂*(I₁+I₂+Iᵤ) + P(Θ₂)*Sᵤ - λ₂*S₂
    du[5] = dI₂ = β₂*S₂*(I₁+I₂+Iᵤ) - (γ+λ₂)*I₂ #+ P(Θ₂)*Iᵤ
    du[6] = dR₂ = γ*I₂ - ρ*R₂ + P(Θ₂)*Rᵤ - λ₂*R₂
    du[7] = dSᵤ = ρ*Rᵤ - βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - (P(Θ₁)+P(Θ₂))*Sᵤ + λ₁*S₁+λ₂*S₂
    du[8] = dIᵤ = βᵤ*Sᵤ*(I₁+I₂+Iᵤ) - γ*Iᵤ + λ₁*I₁ + λ₂*I₂ #- (P(Θ₁)+P(Θ₂))*Iᵤ
    du[9] = dRᵤ = γ*Iᵤ - ρ*Rᵤ - (P(Θ₁)+P(Θ₂))*Rᵤ + λ₁*R₁ + λ₂*R₂
end
tspan = (0,T)

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]

# ρ sweep
for n=1:11
    ρ = 0.002*n
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
    prho = [0.5, 0.4, 0.3, 0.2, 1/7, 0.0033, 0.0033, 0.072, ρ, 0.71, 0.1, 0.1, 1.0, 0.02, 1.1, 1.0, 0.0]
    p=prho
    function save_theta!(u, t, integrator)
        α,β₁,β₂,βᵤ,γ,λ,μ,ρ,b,c₁,c₂,d,r,y₁,y₂,bg = p
        S₁, I₁, R₁, S₂, I₂, R₂, Sᵤ, Iᵤ, Rᵤ = u
    
        I = I₁ + I₂ + Iᵤ
        τ = Sᵤ/(Sᵤ+Rᵤ)
    
        # Define the system of equations
        function equations(F, vars)
    
            Ves, Vgs, Vei, Vgi, Ver, Vgr, Vus, Vui, Vur, Θ₁, Θ₂ = vars
    
            # Wages
            Vfes = (Ves-Vus)*(1 − α)/α
            Vfgs = (Vgs-Vus)*(1 − α)/α
            Vfei = (Vei-Vui)*(1 − α)/α
            Vfgi = (Vgi-Vui)*(1 − α)/α
            Vfer = (Ver-Vur)*(1 − α)/α
            Vfgr = (Vgr-Vur)*(1 − α)/α
            Wes = α*y₁ + (1 - α)*b + (1 − α)*(β₁ - βᵤ)*I*(Vus - Vui) + α*(P(Θ₁)*Vfes + P(Θ₂)*Vfgs)
            Wer = α*y₁ + (1 - α)*b + α*(P(Θ₁)*Vfer + P(Θ₂)*Vfgr)
            Wgs = α*y₂ + (1 - α)*b + (1 − α)*(β₂ - βᵤ)*I*(Vus - Vui) + α*(P(Θ₁)*Vfes + P(Θ₂)*Vfgs)
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
            F[10] = (r + λ)*c₁/Q(Θ₁) + (1 − α)*(b − y₁) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₁ − βᵤ)*I*(Vus − Vui) + (1-τ)*ρ*(Vfer-Vfes) + τ*β₁*I*(Vfes-Vfei)
            F[11] = (r + λ)*c₂/Q(Θ₂) + (1 − α)*(b − y₂) + α*(c₁*Θ₁ + c₂*Θ₂) + τ*(1-α)*(β₂ - βᵤ)*I*(Vus − Vui) + (1-τ)*ρ*(Vfgr-Vfgs) + τ*β₂*I*(Vfgs-Vfgi)
        end
    
        initial_guess = [300, 300, 300, 300, 300, 300, 300, 300, 300, 2.0, 0.01]
        result = nlsolve(equations, initial_guess)
        
        Θ₁ = result.zero[10]
        Θ₂ = result.zero[11]
    
        # Save time and Θ values
        return [Θ₁ Θ₂]
    end
    saved_values = SavedValues(Float64, Array{Float64, 2})
    save_cb = SavingCallback(save_theta!, saved_values, saveat=0:1:T)

    prob = ODEProblem(epiecon_ode, u₀, tspan, prho)
    sol = solve(prob,saveat=0:1:T, callback=save_cb)
    E[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θe[n] = saved_values.saveval[end][1]
    Θg[n] = saved_values.saveval[end][2]
end
plot(0.002:0.002:0.022,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\rho",legend=:outerright)

savefig("sweep_rho.pdf")
