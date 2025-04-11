# McKeithan, T. W. (1995). Kinetic proofreading in T-cell receptor signal transduction. 
#_Proceedings of the national academy of sciences_, _92_(11), 5042-5046.

using SIAN, Logging

# ===================================================================
#               y = TTt)
# ===================================================================

# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff)*C1(t),
    y1(t) = T(t) + C0(t) + C1(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t),
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ R(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff)*C1(t),
    y1(t) = T(t) + C0(t) + C1(t),
    y2(t) = C1(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t),
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t),
    y2(t) = C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff)*C1(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = C2(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kon ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff)*C1(t),
    kon'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t),
    y2(t) = kon(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff)*C1(t),
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t),
    y2(t) = kp(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t),
    kp'(t) = 0,
    y1(t) = T(t)+ C0(t) + C1(t) + C2(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________  EC50 ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t),
    y1(t) = koff/kon + ((T(t)+C0(t)+C1(t)+C2(t)+C3(t))/2),
    y2(t) = T(t)+C0(t)+C1(t)+C2(t)+C3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff)*C3(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((T(t)+C0(t)+C1(t)+C2(t)+C3(t))/2),
    y2(t) = T(t)+C0(t)+C1(t)+C2(t)+C3(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
