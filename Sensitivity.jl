using GlobalSensitivity, Statistics, QuasiMonteCarlo
using Plots, LaTeXStrings
theme(:default,
    lw = 2,
    size = (450,200),
    fontfamily = "computer modern",)

include("ExpandedModel.jl")


T = 1000.0

# Initial Conditions
u₀ = [0.905, 0.0, 0.0, 0.07, 0.0, 0.0, 0.025, 0.0, 0.0]
# Economic Variable ICs [Θ₁,Θ₂,Vfes,Vfgs,Vfei,Vfgi,Vfer,Vfgi]
v₀ = [2.5, 0.001, 25, 25, 25, 25, 25, 25] # Initial guess

v = nlsolve((F,x) -> equations!(F,x,p,u₀),v₀,ftol=1e-14)
u = [u₀;v.zero]
prob = ODEProblem(DAE, u, (0.0,T), p)
#solve to steady state
# Introduce an infection
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0]

# Parameters to change
keys = (:b,:λ)
lbs = [0.01,0.0001] # Lower Bounds
ubs = [0.75,0.01] # Upper Bounds


f1 = function(pars)
    p_new = merge(p,NamedTuple{keys}(pars))
    v = nlsolve((F,x) -> equations!(F,x,p_new,u₀),v₀,ftol=1e-14)
    u = [u₀;v.zero]
    prob = ODEProblem(DAE, u,(0.0,T),p_new)
    sol = solve(prob,saveat = 0:T)
    U = sum(sol[[7,8,9],:],dims=1)
    return [maximum(U),U[end]]
end


A,B = QuasiMonteCarlo.generate_design_matrices(1000,lbs,ubs,LatinHypercubeSample())

SobolIndices = gsa(f1,Sobol(order=[0,1,2]),A,B)

SobolIndices.S2
b1 = bar([L"b",L"\lambda"],SobolIndices.S2[:], title = L"Effect on $U_{max}$", ylims=(0.0,1.0),legend=nothing)
b2 = bar([L"b",L"\lambda"],SobolIndices.S2[2,:], title = L"Effect on $U_{\infty}$",ylims=(0.0,1.4),legend = nothing)
plot(b1,b2, layout = (1,2))




