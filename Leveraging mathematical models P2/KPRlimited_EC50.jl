# Lever, M., Maini, P. K., Van Der Merwe, P. A., & Dushek, O. (2014). 
# Phenotypic models of T cell activation. Nature Reviews Immunology, _14_(9), 619-629.


using SIAN, Logging

#   ANÁLISIS POR SALIDA:
#       1. y = EC50
#           a) Sin conocer params.
#           b) Conociendo kon, koff, T(t)
#           c) Conociendo T(t)
#           d) Conociendo kon
#           e) Conociendo kp
#       2. y = Emax
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+T(t))/2)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))/2)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+T(t))/2)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)+T(t))/2)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ Conociendo kon, kp, T(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+T(t))/2),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

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
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))/2),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t),
    C4'(t) = kp(t)*C3(t) - (koff + phi)*C4(t),
    C5'(t) = phi*C4(t) - (koff)*C5(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+T(t))/2),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t),
    C5'(t) = kp(t)*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - (koff)*C6(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)+T(t))/2),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Conociendo T(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+T(t))/2),
    y2(t) = T(t)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))/2),  
    y2(t) = T(t)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = T(t)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+T(t))/2), 
    y2(t) = T(t)
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)+T(t))/2),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+T(t))/2),
    y2(t) = kon(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))/2),  
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = kon(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + phi)*C4(t),
    C5'(t) = phi*C4(t) - (koff)*C5(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+T(t))/2), 
    y2(t) = kon(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t),
    C5'(t) = kp*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - (koff)*C6(t),
    kon'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)+T(t))/2),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kp'(t) = 0,
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+T(t))/2),
    y2(t) = kp(t)
)

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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))/2),  
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kp'(t) = 0,
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t),
    C4'(t) = kp(t)*C3(t) - (koff + phi)*C4(t),
    C5'(t) = phi*C4(t) - (koff)*C5(t),
    kp'(t) = 0,
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+T(t))/2), 
    y2(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t) +(koff)*C5(t) +(koff)*C6(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t),
    C5'(t) = kp(t)*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - (koff)*C6(t),
    kp'(t) = 0,
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)+T(t))/2),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ T_T ________________________________________________________
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
    y1(t) = koff/kon + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T_T y todos ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = koff/kon(t) + ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))/2),
    y2(t) = C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y = Emax
# ===================================================================

# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff)) * ((C0(t)+C1(t)+C2(t)+T(t)))
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
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^2 * ((C0(t)+C1(t)+C2(t)+C3(t)+T(t)))
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
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t)))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ Conociendo kon, kp, T(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff)) * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

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
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff))^2 * ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Conociendo T(t) ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff)) * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = T(t)
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
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^2 * ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))),
    y2(t) = T(t)
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
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kon'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff)) * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = kon(t)
)

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + phi)*C2(t),
    C3'(t) = phi*C2(t) - (koff)*C3(t),
    kon'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^2 * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^3 * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi)*C1(t),
    C2'(t) = phi*C1(t) - (koff)*C2(t),
    kp'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff)) * ((C0(t)+C1(t)+C2(t)+T(t))),
    y2(t) = kp(t)
)

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
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff))^2 * ((C0(t)+C1(t)+C2(t)+C3(t)+T(t))),
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kp'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ T_T ________________________________________________________
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
    y1(t) = (koff / (koff + phi)) * (kp / (kp + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))),
    y2(t) = C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T_T y todos ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) +(koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + phi)*C3(t),
    C4'(t) = phi*C3(t) - (koff)*C4(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = (koff / (koff + phi)) * (kp(t) / (kp(t) + koff))^3 * ((C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t))),
    y2(t) = C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
