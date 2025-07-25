## Vary parameters

using Plots, LaTeXStrings

using Plots, LaTeXStrings
theme(:wong, lw = 2,
     size = (600,320),
     fontfamily = "Computer Modern",
     )

=======
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

<<<<<<< HEAD
bs = 0.05:0.01:0.75
bsweep = sweep("b",bs)
sweep_b=plot(bs,bsweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"b",legend=:outerright,yaxis=:log2)

ds = 0.01:0.01:2
dsweep = sweep("d",ds)
sweep_d=plot(ds,dsweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"d",legend=:outerright,yaxis=:log2)

β2s = 0.2:0.001:0.4
βg_sweep = sweep("β₂",β2s)
sweep_betag=plot(β2s,βg_sweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_g",legend=:outerright,yaxis=:log2)

β1s = 0.3:0.001:0.6
βe_sweep = sweep("β₁",β1s)
sweep_betae= plot(β1s,βe_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_e",legend=:outerright,yaxis=:log2)

βus = 0.1:0.001:0.3
βu_sweep = sweep("βᵤ",βus)
sweep_betau= plot(βus,βu_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_u",legend=:outerright,yaxis=:log2)

γs = 0.01:0.001:0.25
γ_sweep = sweep("γ",γs)
sweep_gamma = plot(γs,γ_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\gamma",legend=:outerright,yaxis=:log2)

ρs = 0.002:0.00001:0.025
ρ_sweep = sweep("ρ",ρs)
sweep_rho = plot(ρs,ρ_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\rho",legend=:outerright,yaxis=:log2)

plot(sweep_b,sweep_d,sweep_betae,sweep_betag,sweep_betau,sweep_gamma,sweep_rho,layout=(4,2), size = (600,600))
=======
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
>>>>>>> parent of dd1b588 (Sweeps)
savefig("sweep.pdf")
