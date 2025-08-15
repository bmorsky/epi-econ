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
Wf_DF = W₁(0,Θf_DF,Θg_DF,p)
Wg_DF = W₂(0,Θf_DF,Θg_DF,p)

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
Wf = W₁(sol[2,:] ./ F,Θf,Θg,p)
Wg = W₂(sol[5,:] ./ G,Θf,Θg,p)
t = sol.t

pf=plot(t,[F_DF,F],label=[L"F^{DF}" L"F"],xlabel="",legend=:outerright)
pg=plot(t,[G_DF,G],label=[L"G^{DF}" L"G"],xlabel="",legend=:outerright)
pu=plot(t,[U_DF,U],label=[L"U^{DF}" L"U"],xlabel="",legend=:outerright)
psir=plot(t,[S,I,R],label=[L"S" L"I" L"R"],xlabel="",legend=:outerright)
ptightE=plot(t,[Θf_DF,Θf],label=[L"\Theta_f^{DF}" L"\Theta_f"],xlabel="",legend=:outerright)
ptightG=plot(t,[Θg_DF,Θg],label=[L"\Theta_g^{DF}" L"\Theta_g"],xlabel="",legend=:outerright)
pwf=plot(t,[Wf_DF,Wf],label=[L"W_f^{DF}" L"W_f"],xlabel="time (days)",legend=:outerright)
pwg=plot(t,[Wg_DF,Wg],label=[L"W_g^{DF}" L"W_g"],xlabel="time (days)",legend=:outerright)

plot(pf,pg,pu,psir,ptightE,ptightG,pwf,pwg,layout=(4,2))
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
Wf_gig = W₁(sol[2,:] ./ F_gig,Θf_gig,Θg_gig,p_gig)
Wg_gig = W₂(sol[5,:] ./ G_gig,Θf_gig,Θg_gig,p_gig)
t = sol.t

pf_gig=plot(t,[F_DF,F,F_gig],label=[L"F^{DF}" L"F" L"F^B"],xlabel="",legend=:outerright)
pg_gig=plot(t,[G_DF,G,G_gig],label=[L"G^{DF}" L"G" L"G^B"],xlabel="",legend=:outerright)
pu_gig=plot(t,[U_DF,U,U_gig],label=[L"U^{DF}" L"U" L"U^B"],xlabel="",legend=:outerright)
psir_gig=plot(t,[S_gig,I_gig,R_gig],label=[L"\% \Delta S" L"\% \Delta I" L"\% \Delta R"],xlabel="",legend=:outerright)
ptightF_gig=plot(t,[Θf_DF,Θf,Θf_gig],label=[L"\Theta_f^{DF}" L"\Theta_f" L"\Theta_f^B"],xlabel="time (days)",legend=:outerright)
ptightG_gig=plot(t,[Θg_DF,Θg,Θg_gig],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^B"],xlabel="time (days)",legend=:outerright)
pwf_gig=plot(t,[Wf_DF,Wf_gig],label=[L"W_f^{DF}" L"W_f^B"],xlabel="",legend=:outerright)
pwg_gig=plot(t,[Wg_DF,Wg_gig],label=[L"W_g^{DF}" L"W_g^B"],xlabel="",legend=:outerright)

plot(pf_gig,pg_gig,pu_gig,psir_gig,ptightF_gig,ptightG_gig,pwf_gig,pwg_gig,layout=(4,2))
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
Wf_jloss = W₁(sol[2,:] ./ F_jloss,Θf_jloss,Θg_jloss,p_jloss)
Wg_jloss = W₂(sol[5,:] ./ G_jloss,Θf_jloss,Θg_jloss,p_jloss)
t = sol.t

pf_jloss=plot(t,[F_DF,F,F_jloss],label=[L"F^{DF}" L"F" L"F^L"],xlabel="",legend=:outerright,)
pg_jloss=plot(t,[G_DF,G,G_jloss],label=[L"G^{DF}" L"G" L"G^L"],xlabel="",legend=:outerright)
pu_jloss=plot(t,[U_DF,U,U_jloss],label=[L"U^{DF}" L"U" L"U^L"],xlabel="",legend=:outerright)
psir_jloss=plot(t,[S_jloss,I_jloss,R_jloss],label=[L"\% \Delta S" L"\% \Delta I" L"\% \Delta R"],xlabel="",legend=:outerright)
ptightF_jloss=plot(t,[Θf_DF,Θf,Θf_jloss],label=[L"\Theta_f^{DF}" L"\Theta_f" L"\Theta_f^L"],xlabel="time (days)",legend=:outerright)
ptightG_jloss=plot(t,[Θg_DF,Θg,Θg_jloss],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^L"],xlabel="time (days)",legend=:outerright)
pwf_jloss=plot(t,[Wf_DF,Wf_jloss],label=[L"W_f^{DF}" L"W_f^L"],xlabel="",legend=:outerright)
pwg_jloss=plot(t,[Wg_DF,Wg_jloss],label=[L"W_g^{DF}" L"W_g^L"],xlabel="",legend=:outerright)

plot(pf_jloss,pg_jloss,pu_jloss,psir_jloss,ptightF_jloss,ptightG_jloss,pwf_jloss,pwg_jloss,layout=(4,2))
savefig("timeseries_lambdaf.pdf")

# Plot the effect of a higher job losses in the formal economy
prob = ODEProblem(epiecon, u₀, tspan, p_bu)
sol = solve(prob,saveat=0.1, RadauIIA5(), reltol=1e-8, abstol=1e-8)
F_bu = sol[1,:]+sol[2,:]+sol[3,:]
G_bu = sol[4,:]+sol[5,:]+sol[6,:]
U_bu = 1.0 .- F_bu .- G_bu
Θf_bu = sol[10,:]
Θg_bu = sol[11,:]
S_bu = sol[1,:]+sol[4,:]+sol[7,:] - S
I_bu = sol[2,:]+sol[5,:]+sol[8,:] - I
R_bu = - S_bu - I_bu
S_bu = 100*S_bu ./ S
I_bu = 100*I_bu ./ I
R_bu = 100*R_bu ./ R
Wf_bu = W₁(sol[2,:] ./ F_bu,Θf_bu,Θg_bu,p_bu)
Wg_bu = W₂(sol[5,:] ./ G_bu,Θf_bu,Θg_bu,p_bu)
t = sol.t

pf_bu=plot(t,[F_DF,F,F_bu],label=[L"F^{DF}" L"F" L"F^B"],xlabel="",legend=:outerright,)
pg_bu=plot(t,[G_DF,G,G_bu],label=[L"G^{DF}" L"G" L"G^B"],xlabel="",legend=:outerright)
pu_bu=plot(t,[U_DF,U,U_bu],label=[L"U^{DF}" L"U" L"U^B"],xlabel="",legend=:outerright)
psir_bu=plot(t,[S_bu,I_bu,R_bu],label=[L"\% \Delta S" L"\% \Delta I" L"\% \Delta R"],xlabel="",legend=:outerright)
ptightF_bu=plot(t,[Θf_DF,Θf,Θf_bu],label=[L"\Theta_f^{DF}" L"\Theta_f" L"\Theta_f^B"],xlabel="time (days)",legend=:outerright)
ptightG_bu=plot(t,[Θg_DF,Θg,Θg_bu],label=[L"\Theta_g^{DF}" L"\Theta_g" L"\Theta_g^B"],xlabel="time (days)",legend=:outerright)
pwf_bu=plot(t,[Wf_DF,Wf_bu],label=[L"W_f^{DF}" L"W_f^B"],xlabel="",legend=:outerright)
pwg_bu=plot(t,[Wg_DF,Wg_bu],label=[L"W_g^{DF}" L"W_g^B"],xlabel="",legend=:outerright)

plot(pf_bu,pg_bu,pu_bu,psir_bu,ptightF_bu,ptightG_bu,pwf_bu,pwg_bu,layout=(4,2))
savefig("timeseries_bu.pdf")