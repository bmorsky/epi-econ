## Vary parameters

using Plots, LaTeXStrings

<<<<<<< HEAD
<<<<<<< HEAD
using Plots, LaTeXStrings
theme(:wong, lw = 2,
     size = (600,320),
     fontfamily = "Computer Modern",
     )

=======
>>>>>>> parent of dd1b588 (Sweeps)
=======
>>>>>>> parent of dd1b588 (Sweeps)
# Output
E = zeros(101)
G = zeros(101)
U = zeros(101)
Θe = zeros(101)
Θg = zeros(101)

# Parameters
include("model.jl")

T = 2000 # final time
tspan = (0,T)

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]

# b sweep
for n=1:101
    b = 0.75*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b, c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pb = merge(p,(;b=b))

    prob = ODEProblem(epiecon, u₀, tspan, pb)
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
    pbetag= merge(p,(;β₂=β₂))

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
    pd = merge(p,(;d=d))

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
