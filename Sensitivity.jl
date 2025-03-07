using GlobalSensitivity, Statistics, QuasiMonteCarlo
using Plots, LaTeXStrings


include("model.jl")


T = 1000.0

# 1% initially infected
u₀ = [0.905*0.99, 0.905*0.01, 0.0, 0.07*0.99, 0.07*0.01, 0.0, 0.025*0.99, 0.025*0.01, 0.0, 2.8, 0.0165]

prob = ODEProblem(epiecon, u₀,(0.0,T), p)
#solve to steady state

# Parameters to change
keys = ["b","d","bg","λ₁","λ₂"]
lbs = [0.1,0.25,0.0,0.0033,0.0033] # Lower Bounds
ubs = [0.75,1.0,0.09,0.01,0.1] # Upper Bounds


f1 = function(pars)
    p_new = merge(p,NamedTuple{(Symbol.(keys)...,)}(pars))
    prob = ODEProblem(epiecon, u₀,(0.0,T),p_new)
    sol = solve(prob,saveat=1,reltol=1e-8)
    U = sum(sol[[7,8,9],:],dims=1)
    Is = sum(sol[[2,5,8],:],dims=1)
    return [maximum(U),U[end],maximum(Is),Is[end]]
end


A,B = QuasiMonteCarlo.generate_design_matrices(1000,lbs,ubs,LatinHypercubeSample())

SobolIndices = gsa(f1,Sobol(order=[0,1]),A,B)

indices = [L"b",L"d",L"bg",L"\lambda_1",L"$\lambda_2$"]

b1 = bar(indices,SobolIndices.S1[1,:], ylabel = "First order", title = L"$U_{max}$ Sobol Indices",legend=nothing, ylims = (0.0,1.0));
b2 = bar(indices,SobolIndices.S1[2,:], title = L"$U_{\infty}$ Sobol Indices",legend = nothing,ylims=(0.0,1.0));
b3 = bar(indices,SobolIndices.ST[1,:], ylabel = "Total order", legend=nothing,ylims=(0.0,1.0));
b4 = bar(indices,SobolIndices.ST[2,:], legend=nothing,ylims=(0.0,1.0));
plot(b1,b2,b3,b4 ,layout = (2,2))

c1 = bar(indices,SobolIndices.S1[3,:], ylabel = "First order", title = L"$I_{max}$ Sobol Indices",legend=nothing, ylims = (0.0,1.0));
c2 = bar(indices,SobolIndices.S1[4,:], title = L"$I_{\infty}$ Sobol Indices",legend = nothing,ylims=(0.0,1.0));
c3 = bar(indices,SobolIndices.ST[3,:], ylabel = "Total order", legend=nothing,ylims=(0.0,1.0));
c4 = bar(indices,SobolIndices.ST[4,:], legend=nothing,ylims=(0.0,1.0));
plot(c1,c2,c3,c4, layout = (2,2))


