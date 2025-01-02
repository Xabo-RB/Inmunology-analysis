#François, P., Voisinne, G., Siggia, E. D., Altan-Bonnet, G., & Vergassola, M. (2013). 
#Phenotypic model for early T-cell activation displaying sensitivity, specificity, and antagonism. 
#Proceedings of the National Academy of Sciences, 110(10), E888-E897.

using SIAN, Logging

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
    y1(t) = T(t)  + C0(t) + C1(t) + C2(t)

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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t)
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
    y1(t) = T(t) + C0(t)+ C1(t) + C2(t) + C3(t) + C4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ R(t) ________________________________________________________

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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + C4(t),
    y2(t) = C4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + C4(t),
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ kp ________________________________________________________
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ conociendo TODOS ________________________________________________________
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t),
    y2(t) = C3(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y =  Emax 
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
    y1(t) = (kp * (koff^2 + koff * kp + kp^2 + (b) * (koff + kp))) / (kp * (koff + kp)^2 + (b)^2 * (koff + 2 * kp) + (b) * (koff^2 + 4 * koff * kp + 2 * kp^2)) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t)________________________________________________________
# -------------------  N = 2 -------------------

ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t),
    y1(t) = (kp * (koff^2 + koff * kp + kp^2 + (b + gamma * S(t)) * (koff + kp))) / (kp * (koff + kp)^2 + (b + gamma * S(t))^2 * (koff + 2 * kp) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp + 2 * kp^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
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
    y1(t) = (kp * (koff^2 + koff * kp + kp^2 + (b + gamma * S(t)) * (koff + kp))) / (kp * (koff + kp)^2 + (b + gamma * S(t))^2 * (koff + 2 * kp) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp + 2 * kp^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp________________________________________________________
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
    y1(t) = (kp(t) * (koff^2 + koff * kp(t) + kp(t)^2 + (b + gamma * S(t)) * (koff + kp(t)))) / (kp(t) * (koff + kp(t))^2 + (b + gamma * S(t))^2 * (koff + 2 * kp(t)) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp(t) + 2 * kp(t)^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp T(t)________________________________________________________
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
    y1(t) = (kp(t) * (koff^2 + koff * kp(t) + kp(t)^2 + (b + gamma * S(t)) * (koff + kp(t)))) / (kp(t) * (koff + kp(t))^2 + (b + gamma * S(t))^2 * (koff + 2 * kp(t)) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp(t) + 2 * kp(t)^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kon(t),
    y3(t) = T(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Tt________________________________________________________
# -------------------  N = 2 -------------------

ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t),
    y1(t) = (kp * (koff^2 + koff * kp + kp^2 + (b + gamma * S(t)) * (koff + kp))) / (kp * (koff + kp)^2 + (b + gamma * S(t))^2 * (koff + 2 * kp) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp + 2 * kp^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp Tt________________________________________________________
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
    y1(t) = (kp * (koff^2 + koff * kp + kp^2 + (b + gamma * S(t)) * (koff + kp))) / (kp * (koff + kp)^2 + (b + gamma * S(t))^2 * (koff + 2 * kp) + (b + gamma * S(t)) * (koff^2 + 4 * koff * kp + 2 * kp^2)) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y =  EC50 
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
    y1(t) =  (kon/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t)________________________________________________________
# -------------------  N = 2 -------------------

ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t),
    y1(t) =  (kon/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
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
    y1(t) =  (kon(t)/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp________________________________________________________
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
    y1(t) =  (kon/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon T(t) kp________________________________________________________
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
    y1(t) =  (kon(t)/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
    y2(t) = kon(t),
    y3(t) = T(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Tt________________________________________________________
# -------------------  N = 2 -------------------

ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = kp*C0(t) - (koff + kp + b + gamma*S(t))*C1(t) + (b + gamma*S(t))*C2(t),
    C2'(t) = kp*C1(t) - (koff + b + gamma*S(t))*C2(t),
    S'(t) = alpha*C1(t)*(ST - S(t)) - beta*S(t),
    y1(t) =  (kon/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
    y2(t) = T(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp Tt________________________________________________________
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
    y1(t) =  (kon(t)/koff) + (T(t) + C0(t) + C1(t) + C2(t))/2,
    y2(t) = T(t) + C0(t) + C1(t) + C2(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))