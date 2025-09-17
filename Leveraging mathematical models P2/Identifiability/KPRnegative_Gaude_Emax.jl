# Gauud, G., Achar, S., Bourassa, F. X., Davies, J., Hatzihristidis, T., Choi, S., ... & Love, P. E. (2023). 
#CD3ζ ITAMs enable ligand discrimination and antagonism by inhibiting TCR signaling in response to low-affinity 
#peptides. _Nature Immunology_, _24_(12), 2121-2134.

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
    y1(t) = (k2 * (b + koff + kp) * (b^2 + 2*b*koff + koff*kp) + (kp * (b + koff) + koff * (b + koff) + kp^2) * (k1 * (b + koff) + kp))/(b^2 * (koff + 2*kp) + b*(koff^2 + 4*koff*kp + 2*kp^2) + kp * (koff + kp)^2)*(T(t) + C0(t) + C1(t) + C2(t))
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
    y1(t) = (k2 * (b + koff + kp) * (b^2 + 2*b*koff + koff*kp) + (kp * (b + koff) + koff * (b + koff) + kp^2) * (k1 * (b + koff) + kp))/(b^2 * (koff + 2*kp) + b*(koff^2 + 4*koff*kp + 2*kp^2) + kp * (koff + kp)^2)*(T(t) + C0(t) + C1(t) + C2(t)),
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
    y1(t) =  (k2 * (b + koff + kp) * (b^2 + 2*b*koff + koff*kp) + (kp * (b + koff) + koff * (b + koff) + kp^2) * (k1 * (b + koff) + kp))/(b^2 * (koff + 2*kp) + b*(koff^2 + 4*koff*kp + 2*kp^2) + kp * (koff + kp)^2)*(T(t) + C0(t) + C1(t) + C2(t)),
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
    y1(t) =  (k2 * (b + koff + kp) * (b^2 + 2*b*koff + koff*kp) + (kp * (b + koff) + koff * (b + koff) + kp^2) * (k1 * (b + koff) + kp))/(b^2 * (koff + 2*kp) + b*(koff^2 + 4*koff*kp + 2*kp^2) + kp * (koff + kp)^2)*(T(t) + C0(t) + C1(t) + C2(t)),
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
    y1(t) = (k2 * (b + koff + kp(t)) * (b^2 + 2*b*koff + koff*kp(t)) + (kp(t) * (b + koff) + koff * (b + koff) + kp(t)^2) * (k1 * (b + koff) + kp(t)))/(b^2 * (koff + 2*kp(t)) + b*(koff^2 + 4*koff*kp(t) + 2*kp(t)^2) + kp(t) * (koff + kp(t))^2)*(T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y =  T_T
# ===================================================================

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
    y1(t) = T(t)  + C0(t) + C1(t) + C2(t),
    y2(t) =  (k2 * (b + koff + kp) * (b^2 + 2*b*koff + koff*kp) + (kp * (b + koff) + koff * (b + koff) + kp^2) * (k1 * (b + koff) + kp))/(b^2 * (koff + 2*kp) + b*(koff^2 + 4*koff*kp + 2*kp^2) + kp * (koff + kp)^2)*(T(t) + C0(t) + C1(t) + C2(t))
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
    y1(t) = T(t)  + C0(t) + C1(t) + C2(t),
    y2(t) = (k2 * (b + koff + kp(t)) * (b^2 + 2*b*koff + koff*kp(t)) + (kp(t) * (b + koff) + koff * (b + koff) + kp(t)^2) * (k1 * (b + koff) + kp(t)))/(b^2 * (koff + 2*kp(t)) + b*(koff^2 + 4*koff*kp(t) + 2*kp(t)^2) + kp(t) * (koff + kp(t))^2)*(T(t) + C0(t) + C1(t) + C2(t)),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
