# Andersen, P. S., Menné, C., Mariuzza, R. A., Geisler, C., & Karjalainen, K. (2001). 
#A response calculus for immobilized T cell receptor ligands. _Journal of Biological Chemistry_, _276_(52), 49125-49132.

using SIAN, Logging

# ============================================
#               y = C0(t) = R(t) 
# ============================================

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = C0(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = C0(t),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = C0(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = C0(t),
    y2(t) = T(t),
    y3(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y = koff/kon + 0.5*T(t)*C(t) = EC50 
# ===================================================================

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = koff/kon + (T(t)*C0(t))/2
)

ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + (T(t)*C0(t))/2,
    y2(t) = T(t),
    y3(t) = kon(t)
)

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = koff/kon + (T(t)*C0(t))/2,
    y2(t) = T(t)
)

ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + (T(t)*C0(t))/2,
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ============================================
#               y = T(t) 
# ============================================

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = T(t)
)

ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ============================================
#               y = Tt(t) 
# ============================================

ode = @ODEmodel(
    C0'(t) = kon * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t),
    y1(t) = T(t) + C0(t)
)

ode = @ODEmodel(
    C0'(t) = kon(t) * P(t) * T(t) - koff*C0(t),
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t),
    kon'(t) = 0,
    y1(t) = T(t) + C0(t),
    y3(t) = kon(t)
)



@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

 