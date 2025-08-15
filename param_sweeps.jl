## Vary parameters

using Plots, LaTeXStrings
theme(:wong, lw = 2,
     size = (600,320),
     fontfamily = "Computer Modern",
     )

# Output
F = zeros(101)
G = zeros(101)
U = zeros(101)
Θf = zeros(101)
Θg = zeros(101)
Wf = zeros(101)
Wg = zeros(101)

# Parameters
include("model.jl")

T = 2000 # final time
tspan = (0,T)

# 1% initially infected
u₀ = [SF₀*0.99, SF₀*0.01, 0.0, SG₀*0.99, SG₀*0.01, 0.0, SU₀*0.99, SU₀*0.01, 0.0, Θ₁₀, Θ₂₀]

# bᵤ sweep
for n=1:101
    bu = 0.35+0.1*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b, c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pbu = merge(p,(;bᵤ=bu))

    prob = ODEProblem(epiecon, u₀, tspan, pbu)
    sol = solve(prob,saveat=1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
    F[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θf[n] = sol[10,end]
    Θg[n] = sol[11,end]
    Wf[n] = W₁(sol[2,end] / F[n],Θf[n],Θg[n],pbu)
    Wg[n] = W₂(sol[5,end] / G[n],Θf[n],Θg[n],pbu)
end
sweep_bu=plot(collect(LinRange(0.35,0.45,101)),[F,G,U,Θf,Wf,Wg],label=[L"F" L"G" L"U" L"\Theta_f" L"W_f" L"W_g"],xlabel=L"b_u",legend=:outerright)

# b₂ sweep
for n=1:101
    bg = 0.15+0.1*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b, c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pbu = merge(p,(;b₂=bg))

    prob = ODEProblem(epiecon, u₀, tspan, pbu)
    sol = solve(prob,saveat=1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
    F[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θf[n] = sol[10,end]
    Θg[n] = sol[11,end]
    Wf[n] = W₁(sol[2,end] / F[n],Θf[n],Θg[n],pbu)
    Wg[n] = W₂(sol[5,end] / G[n],Θf[n],Θg[n],pbu)
end
sweep_bg=plot(collect(LinRange(0.15,0.25,101)),[F,G,U,Θf,Wf,Wg],label=[L"F" L"G" L"U" L"\Theta_f" L"W_f" L"W_g"],xlabel=L"b_g",legend=:outerright)

# β₂ sweep
for n=1:101
    βg = 0.34+0.6*(n-1)/100
    #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
    pbetag= merge(p,(;β₂=βg))

    prob = ODEProblem(epiecon, u₀, tspan, pbetag)
    sol = solve(prob,saveat=1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
    F[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θf[n] = sol[10,end]
    Θg[n] = sol[11,end]
    Wf[n] = W₁(sol[2,end] / F[n],Θf[n],Θg[n],pbetag)
    Wg[n] = W₂(sol[5,end] / G[n],Θf[n],Θg[n],pbetag)
end
sweep_betag=plot(collect(LinRange(0.34,0.4,101)),[F,G,U,Θf,Wf,Wg],label=[L"F" L"G" L"U" L"\Theta_f" L"W_f" L"W_g"],xlabel=L"\beta_g",legend=:outerright)

# # d sweep
# for n=1:101
#     d = (n-1)/100
#     #    α,   β₁,  β₂,  βᵤ,  γ,   λ₁,     λ₂,     μ,     ρ,    b,    c₁,  c₂,  d, r,    y₁,  y₂,  bg
#     pd = merge(p,(;d=d))

#     prob = ODEProblem(epiecon_ode, u₀, tspan, pd)
#     sol = solve(prob,saveat=1)
#     E[n] = sol[1,end]+sol[2,end]+sol[3,end]
#     G[n] = sol[4,end]+sol[5,end]+sol[6,end]
#     U[n] = sol[7,end]+sol[8,end]+sol[9,end]
#     Θe[n] = sol[10,end]
#     Θg[n] = sol[11,end]
# end
# sweep_d=plot(0:0.01:1,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"d",legend=:outerright)

plot(sweep_bu,sweep_bg,sweep_betag,layout=(2,2))
savefig("sweep.pdf")
