## The gig economy during an epidemic: coupling disease transmission with labour market dynamics
#### Bryce Morsky, Tyler Meadows, Felicia M.G. Magpantay, and Troy Day

This repository contains the code used for simulations and figures in "The gig economy during an epidemic: coupling disease transmission with labour market dynamics" (20??). 

<b> Updates (March 7, 2025) </b>


- Cleaned up (removed) some old code and renamed some files. model.jl now contains the main model, epiecon.jl creates the time-series figures for the main model, and epiecon_expanded.jl creates the figure for the expanded model. 
- Fixed some small numerical issues in solving the model. I'm using Julia's Wong theme for plots, since it is colourblind friendly and looks pretty okay. 
- Sensitivity analysis is up and running for the basic model. Comparing the effects of $b,d,bg,\lambda_1,$ and $\lambda_2$ on $U_\infty,U_{max},I_\infty$ and $I_{max}$, where $I_\infty = \lim_{t\to\infty}I(t)$ and $I_{max} = \max \{ I(t): t\ge 0 \} $
- Parameter sweeps organized into a single file and completed for parameters $b,d,\beta_e,\beta_g,\beta_u,\gamma$ and $\rho$.

<b> Updates (Feb 28, 2025) </b>

- Started working on the sensitivity analysis. Ran into some issues with $V_{fgs} and $V_{fgr}$ becoming negative, which caused problems calculating $\Theta_g$. 
- Fixed some typos in the model definition: $V_{us}-V_{ui}$ had $\lambda$ in a few places that should have been $\gamma$.

<b> Updates (Feb 27, 2025) </b>

- Moved the model definition to its own file. You can load this into other files by writing
    ```Julia
    include("ExpandedModel.jl") 
    ```
    This gives access to
    -  the differential algebraic system <code>DAE</code>, which can be solved using the DifferentialEquations.jl ODE solver
    - The algebraic equations <code>equations!(F,vars,p,u)</code> which can be used with <code>nlsolve</code> to find consistent initial conditions for the economic variables
    - A base set of parameters <code>p</code> and the matching functions <code>P(Theta)</code> and <code>Q(Theta)</code>. The parameters can be altered piece by piece by defining new parameters using <code> merge</code>. For example, to change only the value of $b$ in the parameters
    ```Julia 
    new_parameters = merge(p,(;b=0.01))    
    ```
- Fixed parameter sweep file for $b$