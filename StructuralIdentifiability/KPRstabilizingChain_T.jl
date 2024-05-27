#Gálvez, J., Galvez, J. J., & Garcia-Penarrubia, P. (2016). 
#TCR/pMHC interaction: phenotypic model for an unsolved enigma. 
#Frontiers in Immunology, 7, 228066.

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
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = T(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (3/(1 + 2*r))*kp*C2(t) - (4/(1 + 3*r))*koff*C3(t),
    y1(t) = T(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (r^2)*kp*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp)*C3(t),
    C4'(t) = (r^3)*kp*C3(t) - ((5/(1 + 4*r))*koff)*C4(t),
    y1(t) = T(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (r^2)*kp*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp)*C3(t),
    C4'(t) = (r^3)*kp*C3(t) - ((5/(1 + 4*r))*koff + (r^4)*kp)*C4(t),
    C5'(t) = (r^4)*kp*C4(t) - (6/(1 + 5*r))*koff*C5(t),
    y1(t) = T(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (3/(1 + 2*r))*kp(t)*C2(t) - (4/(1 + 3*r))*koff*C3(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (r^2)*kp(t)*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp(t))*C3(t),
    C4'(t) = (r^3)*kp(t)*C3(t) - ((5/(1 + 4*r))*koff)*C4(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (r^2)*kp(t)*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp(t))*C3(t),
    C4'(t) = (r^3)*kp(t)*C3(t) - ((5/(1 + 4*r))*koff + (r^4)*kp(t))*C4(t),
    C5'(t) = (r^4)*kp(t)*C4(t) - (6/(1 + 5*r))*koff*C5(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (3/(1 + 2*r))*kp*C2(t) - (4/(1 + 3*r))*koff*C3(t),
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (r^2)*kp*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp)*C3(t),
    C4'(t) = (r^3)*kp*C3(t) - ((5/(1 + 4*r))*koff)*C4(t),
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + r*kp)*C1(t),
    C2'(t) = r*kp*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp)*C2(t),
    C3'(t) = (r^2)*kp*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp)*C3(t),
    C4'(t) = (r^3)*kp*C3(t) - ((5/(1 + 4*r))*koff + (r^4)*kp)*C4(t),
    C5'(t) = (r^4)*kp*C4(t) - (6/(1 + 5*r))*koff*C5(t),
    y1(t) = T(t),
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (3/(1 + 2*r))*kp(t)*C2(t) - (4/(1 + 3*r))*koff*C3(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (r^2)*kp(t)*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp(t))*C3(t),
    C4'(t) = (r^3)*kp(t)*C3(t) - ((5/(1 + 4*r))*koff)*C4(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t) + ((5/(1 + 4*r))*koff)*C4(t) + (6/(1 + 5*r))*koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + r*kp(t))*C1(t),
    C2'(t) = r*kp(t)*C1(t) - ((3/(1 + 2*r))*koff + (r^2)*kp(t))*C2(t),
    C3'(t) = (r^2)*kp(t)*C2(t) - ((4/(1 + 3*r))*koff + (r^3)*kp(t))*C3(t),
    C4'(t) = (r^3)*kp(t)*C3(t) - ((5/(1 + 4*r))*koff + (r^4)*kp(t))*C4(t),
    C5'(t) = (r^4)*kp(t)*C4(t) - (6/(1 + 5*r))*koff*C5(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

