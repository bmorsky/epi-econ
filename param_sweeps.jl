## Vary parameters

using Plots, LaTeXStrings

# Output

E = Float64[]
G = Float64[]
U = Float64[]
Θe = Float64[]
Θg = Float64[]
# Parameters
include("model.jl")


function sweep(tag,vals)
    T = 2000 # final time
    tspan = (0,T)
    # 1% initially infected
    u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]
    # Initialize storage
    E = Float64[]
    G = Float64[]
    U = Float64[]
    Θe = Float64[]
    Θg = Float64[]
    for val ∈ vals
        pb = merge(p,NamedTuple{(Symbol(tag),)}([val]))
        prob = ODEProblem(epiecon, u₀, tspan, pb)
        sol = solve(prob,saveat=1)
        push!(E,sol[1,end]+sol[2,end]+sol[3,end])
        push!(G,sol[4,end]+sol[5,end]+sol[6,end])
        push!(U,sol[7,end]+sol[8,end]+sol[9,end])
        push!(Θe,sol[10,end])
        push!(Θg,sol[11,end])
    end
    return [E,G,U,Θe,Θg]
end

bs = 0.01:0.01:0.75
bsweep = sweep("b",bs)
sweep_b=plot(bs,bsweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"b",legend=:outerright)

β2s = 0.2:0.001:0.4
βg_sweep = sweep("β₂",β2s)
sweep_betag=plot(β2s,βg_sweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_g",legend=:outerright)

ds = 0.0:0.01:1.0
dsweep= sweep("d",ds)
sweep_d=plot(ds,dsweep,label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"d",legend=:outerright)

β1s = 0.3:0.001:0.6
βe_sweep = sweep("β₁",β1s)
sweep_betae= plot(β1s,βe_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_e",legend=:outerright)

βus = 0.1:0.001:0.3
βu_sweep = sweep("βᵤ",βus)
sweep_betau= plot(βus,βu_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\beta_u",legend=:outerright)

γs = 0.01:0.001:0.25
γ_sweep = sweep("γ",γs)
sweep_gamma = plot(γs,γ_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\gamma",legend=:outerright)

ρs = 0.002:0.00001:0.025
ρ_sweep = sweep("ρ",ρs)
sweep_rho = plot(ρs,ρ_sweep,label = [L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"\rho",legend=:outerright)

plot(sweep_b,sweep_betae,sweep_betag,sweep_betau,sweep_gamma,sweep_rho,layout=(3,2))
savefig("Figures/sweep.pdf")
