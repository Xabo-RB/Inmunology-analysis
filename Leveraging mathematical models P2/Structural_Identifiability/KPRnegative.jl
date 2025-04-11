#François, P., Voisinne, G., Siggia, E. D., Altan-Bonnet, G., & Vergassola, M. (2013). 
#Phenotypic model for early T-cell activation displaying sensitivity, specificity, and antagonism. 
#Proceedings of the National Academy of Sciences, 110(10), E888-E897.

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
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C2(t)
)
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + b + gamma*S(t))*C3(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C3(t)
)
# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + b + gamma*S(t))*C4(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C4(t)
)
# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + kp + b + gamma*S(t))*C4(t) + (b + gamma*S(t))*C5(t),
    C5'(t) = kp*C4(t) - (koff + b + gamma*S(t))*C5(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C5(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
# __________ conociendo TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C2(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + b + gamma*S(t))*C3(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C3(t),

)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t) + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp(t)*C3(t) - (koff + b + gamma*S(t))*C4(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C4(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t) + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t) + b + gamma*S(t))*C4(t) + (b + gamma*S(t))*C5(t),
    C5'(t) = kp(t)*C4(t) - (koff + b + gamma*S(t))*C5(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C5(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C2(t),
    y2(t) = T(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + b + gamma*S(t))*C3(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C3(t),
    y2(t) = T(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + b + gamma*S(t))*C4(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C4(t),
    y2(t) = T(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + kp + b + gamma*S(t))*C4(t) + (b + gamma*S(t))*C5(t),
    C5'(t) = kp*C4(t) - (koff + b + gamma*S(t))*C5(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    y1(t) = C5(t),
    y2(t) = T(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    y1(t) = C2(t),
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + b + gamma*S(t))*C3(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    y1(t) = C3(t),
    y2(t) = kon(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + b + gamma*S(t))*C4(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    y1(t) = C4(t),
    y2(t) = kon(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + kp + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp*C2(t) - (koff + kp + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp*C3(t) - (koff + kp + b + gamma*S(t))*C4(t) + (b + gamma*S(t))*C5(t),
    C5'(t) = kp*C4(t) - (koff + b + gamma*S(t))*C5(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kon'(t) = 0,
    y1(t) = C5(t),
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kp'(t) = 0,
    y1(t) = C2(t),
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + b + gamma*S(t))*C3(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kp'(t) = 0,
    y1(t) = C3(t),
    y2(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t) + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp(t)*C3(t) - (koff + b + gamma*S(t))*C4(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kp'(t) = 0,
    y1(t) = C4(t),
    y2(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t) + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t) + b + gamma*S(t))*C2(t) + (b + gamma*S(t))*C3(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t) + b + gamma*S(t))*C3(t) + (b + gamma*S(t))*C4(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t) + b + gamma*S(t))*C4(t) + (b + gamma*S(t))*C5(t),
    C5'(t) = kp(t)*C4(t) - (koff + b + gamma*S(t))*C5(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t), 
    kp'(t) = 0,
    y1(t) = C5(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
