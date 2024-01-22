# Andersen, P. S., Menné, C., Mariuzza, R. A., Geisler, C., & Karjalainen, K. (2001). 
#A response calculus for immobilized T cell receptor ligands. _Journal of Biological Chemistry_, _276_(52), 49125-49132.

using SIAN, Logging

ode = @ODEmodel(
    C0'(t) = - kon * P * T - koff*C0(t),
    y1(t) = C0(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

ode = @ODEmodel(
    C0'(t) = - kon * P * T(t) - koff*C0(t),
    T'(t) = 0,
    y1(t) = C0(t),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

ode = @ODEmodel(
    C0'(t) = - kon(t) * P * T - koff*C0(t),
    kon'(t) = 0,
    y1(t) = C0(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


ode = @ODEmodel(
    C0'(t) = - kon(t) * P * T(t) - koff*C0(t),
    T'(t) = 0,
    kon'(t) = 0,
    y1(t) = C0(t),
    y2(t) = T(t),
    y3(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))