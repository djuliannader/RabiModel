
using Optim

# -----------------------------------
# Quadratic function with extra parameter
# -----------------------------------
function quadratic(x, y, a)
    return -((x - a)^2 + (y + 1)^2)
end

# -----------------------------------
# Fixed parameter
# -----------------------------------
a = 3.0   # ← this is NOT optimized

# -----------------------------------
# Wrapper for optimizer
# -----------------------------------
objective(v) = -quadratic(v..., a)
# minus sign because Optim minimizes

# -----------------------------------
# Run optimization
# -----------------------------------
initial_guess = [0.0, 0.0]

result = optimize(objective, initial_guess)

println("Best parameters (x,y): ", result.minimizer)
println("Maximum value: ", -result.minimum)
