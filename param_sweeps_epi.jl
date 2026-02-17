## Vary parameters and see how they affect the epidemiological steady states
using Plots, LaTeXStrings

# Parameters
include("model.jl")


function sweep_epi(tag,vals)
    T = 2000 # final time
    tspan = (0,T)
    # 1% initially infected
    u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]
    # Initialize storage
    S = Float64[]
    I = Float64[]
    R = Float64[]
    for val ∈ vals
        pb = merge(p,NamedTuple{(Symbol(tag),)}([val]))
        prob = ODEProblem(epiecon, u₀, tspan, pb)
        sol = solve(prob,saveat=1)
        push!(S,sol[1,end]+sol[4,end]+sol[7,end])
        push!(I,sol[2,end]+sol[5,end]+sol[8,end])
        push!(R,sol[3,end]+sol[6,end]+sol[9,end])
    end
    return [S,I,R]
end

bs = 0.01:0.001:0.75 
sweepb = sweep_epi("b",bs)
sweep_b=plot(bs,sweepb,
            label=[L"S" L"I" L"R"],
            xlabel=L"b",
            legend=:outerright,
            )

αs = 0.43:0.001:0.56
αsweep = sweep_epi("α",αs)
sweep_α = plot(αs,αsweep, 
          label=[L"S" L"I" L"R"],
          xlabel=L"α",
          legend=:outerright,
          )

λs = 0.001:0.00001:0.05
λ1sweep = sweep_epi("λ₁",λs)
sweep_λ1 = plot(λs,λ1sweep, 
            label=[L"S" L"I" L"R"],
            xlabel=L"\lambda_e",
            legend=:outerright,
            )
λ2sweep = sweep_epi("λ₂",λs)
sweep_λ2 = plot(λs,λ2sweep, 
            label=[L"S" L"I" L"R"],
            xlabel=L"\lambda_g",
            legend=:outerright,
            )
ds = 0.0:0.01:1.0
dsweep= sweep_epi("d",ds)
sweep_d=plot(ds,dsweep,
            label=[L"S" L"I" L"R"],
            xlabel=L"d",
            legend=:outerright,
            )
y1s =  1.1:0.01:1.4
y1sweep = sweep_epi("y₁",y1s)



plot(sweep_b,sweep_betae,sweep_betag,sweep_betau,sweep_gamma,sweep_rho,layout=(3,2))
savefig("Figures/sweep_econ.pdf")
savefig("Figures/sweep_econ.png")
