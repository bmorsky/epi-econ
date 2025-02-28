using Plots, LaTeXStrings
using NLsolve, OrdinaryDiffEq
using Parameters
theme(:default,
     size = (450,200),
     fontfamily = "computer modern",
     lw = 2
     )


# Parameters
T = 20000.0 # final time

# Matching parameters and functions
η = 0.5
μ = 0.072
P(Θ) = μ*abs(Θ)^(1-η)
Q(Θ) = μ*abs(Θ)^(-η)

# Import system equations
include("ExpandedModel.jl")

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0]
v₀ = [2.5, 0.001, 25, 25, 25, 25, 25, 25] # Initial guess


# b sweep
b_range = 0.0:0.0075:0.75

E = zeros(length(b_range))
G = similar(E)
U = similar(E)
Θe = similar(E)
Θg = similar(E)


for (n,bs) = enumerate(b_range)
    pb = merge(p, (;b=bs,d=0.01))

    v = nlsolve((F,x) -> equations!(F,x,pb,u₀),v₀,ftol=1e-14)
    u = [u₀;v.zero]

    prob = ODEProblem(DAE, u,(0.0,T), pb)
    sol = solve(prob)
    E[n] = sol[1,end]+sol[2,end]+sol[3,end]
    G[n] = sol[4,end]+sol[5,end]+sol[6,end]
    U[n] = sol[7,end]+sol[8,end]+sol[9,end]
    Θe[n] = sol[10,end]
    Θg[n] = sol[11,end]
end
plot(b_range,[E,G,U,Θe,Θg],label=[L"E" L"G" L"U" L"\Theta_e" L"\Theta_g"],xlabel=L"b",legend=:outerright)

savefig("Figures/sweep_b.pdf")
