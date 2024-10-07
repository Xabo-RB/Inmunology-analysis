using DifferentialEquations
using Plots

include("Modelos.jl")

# Parámetros de la ODE (ajústalos según sea necesario)
p = [10, 1, 0.1, 1, 1, 10]

# Condiciones iniciales (ajústalas según sea necesario)
x0 = [100, 2, 0, 0, 0, 0, 0]

(d, tspan) = (1.0e-11, (0.0,100)); # step size and time interval in days

# Resolver la ODE
prob = ODEProblem(ODEKPC, x0, tspan, p)
sol = solve(prob, dt=d)

# Extraer y graficar el tercer estado (x[3])
t_vals = sol.t
x3_vals = sol[3, :]  # x[3] representa a Tp(t) / Rp(t)

# Graficar
plot(t_vals, x3_vals, label="x[3] (Tp(t) / Rp(t))", xlabel="Tiempo (t)", ylabel="x[3]", lw=2)
title!("Solución del tercer estado de la ODE")