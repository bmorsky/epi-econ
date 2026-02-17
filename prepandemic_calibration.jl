using NonlinearSolve

α = 0.5
λ₁ = 0.001096
λ₂ = 0.005479
μ = 0.01517
b = 0.4
r = 0.000137
η = 0.5
W₁ = 1.0
W₂ = 0.85
b₂ = 0.2

# y₁ = u[1], y₂ = u[2]
# c₁ = u[3]

c₁ = μ*(W₂ + b₂ - W₁)/(λ₂ - λ₁)

f(u, p) =  [-W₁ + α*u[1] + (1-α)*b + α*(c₁ + u[3]*(c₁/u[3])^(1/η)), # Wage eqn 1
            -W₂ - (1-α)*b₂ + α*u[2] + (1-α)*b + α*(c₁ + u[3]*(c₁/u[3])^(1/η)), # Wage eqn 2
            (u[1] - W₁)/(r+λ₁) - (u[2] - W₂)/(r+λ₂) # Free-entry 1
            ]

u₀ = [2,2,0.1]
prob = NonlinearProblem(f, u₀)
y₁, y₂, c₂ = solve(prob)

θ₂ = (c₁/c₂)^(1/η)

# (u[1] - W₁)/(r+λ₁) - u[3]/μ, # Free-entry 1
# (u[2] - W₂)/(r+λ₂) - u[3]/μ # Free-entry 2