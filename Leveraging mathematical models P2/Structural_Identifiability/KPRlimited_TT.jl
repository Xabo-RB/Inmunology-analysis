# Lever, M., Maini, P. K., Van Der Merwe, P. A., & Dushek, O. (2014). 
# Phenotypic models of T cell activation. Nature Reviews Immunology, _14_(9), 619-629.


using SIAN, Logging

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Conociendo R(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    kon'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp _______________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Conociendo todo ________________________________________________________
# -------------------  N = 1 -------------------

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
