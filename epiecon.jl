## Make Figures 2,3,4

using Plots, LaTeXStrings
theme(:wong, lw = 2,
     size = (600,400),
     fontfamily = "Computer Modern",
     )

include("model.jl")

# Disease-free scenario
u₀ = [SF₀, 0.0, 0.0, SG₀, 0.0, 0.0, SU₀, 0.0, 0.0, Θ₁₀, Θ₂₀]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F_DF = sol[1,:]+sol[2,:]+sol[3,:]
G_DF = sol[4,:]+sol[5,:]+sol[6,:]
U_DF = 1.0 .- F_DF .- G_DF
Θf_DF = sol[10,:]
Θg_DF = sol[11,:]

# 1% initially infected
u₀ = [SF₀*0.99, SF₀*0.01, 0.0, SG₀*0.99, SG₀*0.01, 0.0, SU₀*0.99, SU₀*0.01, 0.0, Θ₁₀, Θ₂₀]
prob = ODEProblem(epiecon, u₀, tspan, p)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F = sol[1,:]+sol[2,:]+sol[3,:]
G = sol[4,:]+sol[5,:]+sol[6,:]
U = 1.0 .- F .- G
Θf = sol[10,:]
Θg = sol[11,:]

S = sol[1,:]+sol[4,:]+sol[7,:]
I = sol[2,:]+sol[5,:]+sol[8,:]
R = 1.0 .- S .- I
t = sol.t

pf=plot(t,[F_DF,F],label=[L"F^{DF}" L"F"],xlabel="",legend=:outerright)
pg=plot(t,[G_DF,G],label=[L"G^{DF}" L"G"],xlabel="",legend=:outerright)
pu=plot(t,[U_DF,U],label=[L"U^{DF}" L"U"],xlabel="",legend=:outerright)
psir=plot(t,[S,I,R],label=[L"S" L"I" L"R"],xlabel="",legend=:outerright)
ptightE=plot(t,[Θf_DF,Θf],label=[L"\Theta_f^{DF}" L"\Theta_f"],xlabel="time (days)",legend=:outerright)
ptightG=plot(t,[Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel="time (days)",legend=:outerright)

plot(pf,pg,pu,psir,ptightE,ptightG,layout=(3,2))
savefig("timeseries.pdf")

# Plot the effect of a government benefit to gig workers
prob = ODEProblem(epiecon, u₀, tspan, p_gig)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F_gig = sol[1,:]+sol[2,:]+sol[3,:]
G_gig = sol[4,:]+sol[5,:]+sol[6,:]
U_gig = 1.0 .- F_gig .- G_gig
Θf_gig = sol[10,:]
Θg_gig = sol[11,:]

S_gig = sol[1,:]+sol[4,:]+sol[7,:] - S
I_gig = sol[2,:]+sol[5,:]+sol[8,:] - I
R_gig = - S_gig - I_gig

S_gig = 100*S_gig ./ S
I_gig = 100*I_gig ./ I
R_gig = 100*R_gig ./ R

t = sol.t

pf_gig=plot(t,[F_DF,F,F_gig],label=[L"F^{DF}" L"F" L"F^B"],xlabel="",legend=:outerright)
pg_gig=plot(t,[G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^B"],xlabel="",legend=:outerright)
pu_gig=plot(t,[U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^B"],xlabel="",legend=:outerright)
psir_gig=plot(t,[S_gig,I_gig,R_gig],label=[L"\% \Delta S" L"\% \Delta I" L"\% \Delta R"],xlabel="",legend=:outerright)
ptightF_gig=plot(t,[Θf_DF,Θf,Θf_gig],label=[L"\Theta_f^{DF}" L"\Theta_f" L"\Theta_f^B"],xlabel="time (days)",legend=:outerright)
ptightG_gig=plot(t,[Θg_DF,Θg,Θg_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^B"],xlabel="time (days)",legend=:outerright)

plot(pf_gig,pg_gig,pu_gig,psir_gig,ptightF_gig,ptightG_gig,layout=(3,2))
savefig("timeseries_gig_benefits.pdf")

# Plot the effect of a higher job losses in the formal economy
prob = ODEProblem(epiecon, u₀, tspan, p_jloss)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F_jloss = sol[1,:]+sol[2,:]+sol[3,:]
G_jloss = sol[4,:]+sol[5,:]+sol[6,:]
U_jloss = 1.0 .- F_jloss .- G_jloss
Θf_jloss = sol[10,:]
Θg_jloss = sol[11,:]

S_jloss = sol[1,:]+sol[4,:]+sol[7,:] - S
I_jloss = sol[2,:]+sol[5,:]+sol[8,:] - I
R_jloss = - S_jloss - I_jloss

S_jloss = 100*S_jloss ./ S
I_jloss = 100*I_jloss ./ I
R_jloss = 100*R_jloss ./ R

pf_jloss=plot(t,[F_DF,F,F_jloss],label=[L"F^{DF}" L"F" L"J^L"],xlabel="",legend=:outerright,)
pg_jloss=plot(t,[G_DF,G,G_jloss],label=[L"G^{DF}" L"G" L"G^L"],xlabel="",legend=:outerright)
pu_jloss=plot(t,[U_DF,U,U_jloss],label=[L"U^{DF}" L"U" L"U^L"],xlabel="",legend=:outerright)
psir_jloss=plot(t,[S_jloss,I_jloss,R_jloss],label=[L"\% \Delta S" L"\% \Delta I" L"\% \Delta R"],xlabel="",legend=:outerright)
ptightF_jloss=plot(t,[Θf_DF,Θf,Θf_jloss],label=[L"\Theta_f^{DF}" L"\Theta_f" L"\Theta_f^L"],xlabel="time (days)",legend=:outerright)
ptightG_jloss=plot(t,[Θg_DF,Θg,Θg_jloss],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^L"],xlabel="time (days)",legend=:outerright)

plot(pf_jloss,pg_jloss,pu_jloss,psir_jloss,ptightF_jloss,ptightG_jloss,layout=(3,2))
savefig("timeseries_quarantine.pdf")