## The gig economy during an epidemic: coupling disease transmission with labour market dynamics
#### Bryce Morsky, Tyler Meadows, Felicia M.G. Magpantay, and Troy Day

This repository contains the code used for simulations and figures in "The gig economy during an epidemic: coupling disease transmission with labour market dynamics" (20??). 
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