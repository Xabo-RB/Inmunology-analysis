using DifferentialEquations
using Plots

include("Modelos.jl")

# Parámetros de la ODE
p = complex([10, 1, 0.1, 1, 1, 10])

# Condiciones iniciales
x0 = complex([100, 2, 0, 0, 0, 0, 0])

(d, tspan) = (1.0e-11, (0.0,700)); # step size and time interval in days

# Resolver la ODE
prob = ODEProblem(ODEKPC, x0, tspan, p)
sol = solve(prob, dt=d)


t_vals = sol.t
x3_vals = sol[3, :]

# Graficar
plot(t_vals, x3_vals, label="Response", xlabel="Tiempo (t)", ylabel="x[3]", lw=2)
title!("System Response")