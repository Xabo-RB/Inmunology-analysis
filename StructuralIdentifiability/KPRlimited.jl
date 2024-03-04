# Lever, M., Maini, P. K., Van Der Merwe, P. A., & Dushek, O. (2014). 
# Phenotypic models of T cell activation. Nature Reviews Immunology, _14_(9), 619-629.


using SIAN, Logging

#   ANÁLISIS POR SALIDA:
#       1. y = Cn
#           a) Sin conocer params.
#           b) Conociendo kon, koff, T(t)
#           c) Conociendo T(t)
#           d) Conociendo kon
#           e) Conociendo kp
#       2. y = EC50
#           a) Sin conocer params.
#           b) Conociendo kon, koff, T(t)
#           c) Conociendo T(t)
#           d) Conociendo kon
#           e) Conociendo kp
#       3. y = T(t)
#           a) Sin conocer params.
#           b) Conociendo kon, koff, T(t)
#           c) Conociendo T(t)
#           d) Conociendo kon
#           e) Conociendo kp


# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    y1(t) = C1(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    y1(t) = C2(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    y1(t) = C3(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + phi)*C4(t),
    C5'(t) = phi*C4(t) - (koff)*C5(t),
    y1(t) = C4(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t),
    C5'(t) = kp*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - (koff)*C6(t),
    y1(t) = C5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ Conociendo kon, kp, T(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff(t) + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    y1(t) = C1(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    y1(t) = C2(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    y1(t) = C3(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + phi)*C4(t),
    C5'(t) = phi*C4(t) - (koff)*C5(t),
    y1(t) = C4(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t),
    C5'(t) = kp*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - (koff)*C6(t),
    y1(t) = C5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))