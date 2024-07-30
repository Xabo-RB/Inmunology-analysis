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
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))*D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    y1(t) = C2(t) + D2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ R(t) y T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))*D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    y1(t) = C2(t) + D2(t),
    y2(t) = R - C0(t) - C1(t) - C2(t) - D0(t) - D1(t) - D2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ R(t) y TT ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))*D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    R'(t) = 0,
    y1(t) = C2(t) + D2(t),
    y2(t) = R(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))*D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    y1(t) = R - C0(t) - C1(t) - C2(t) - D0(t) - D1(t) - D2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ conociendo TODOS con Tt________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    k'(t) = 0,
    phi'(t) = 0,
    R'(t) = 0,
    y1(t) = C2(t) + D2(t),
    y2(t) = R(t),
    y3(t) = k(t),
    y4(t) = phi(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo TODOS con T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    k'(t) = 0,
    phi'(t) = 0,
    y1(t) = C2(t) + D2(t),
    y2(t) = R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t),
    y3(t) = k(t),
    y4(t) = phi(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kon y R(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))* C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))* D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    k'(t) = 0,
    y2(t) = C2(t) + D2(t),
    y3(t) = k(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kon y T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))* C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))* D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    k'(t) = 0,
    y2(t) = R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t),
    y3(t) = k(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kp y R(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    phi'(t) = 0,
    y1(t) = C2(t) - D2(t),
    y2(t) = phi(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kp y T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    phi'(t) = 0,
    y1(t) = R - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t),
    y2(t) = phi(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ TT ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))*C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))*D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    R'(t) = 0,
    y1(t) = R(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kp y TT ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    phi'(t) = 0,
    R'(t) = 0,
    y1(t) = R(t),
    y2(t) = phi(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kon y T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi) * C0(t) + (b + gamma*S(t))* C1(t),
    C1'(t) = phi*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi + b + gamma*S(t))*C1(t),
    C2'(t) = phi*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi) * D0(t) + (b + gamma*S(t))* D1(t),
    D1'(t) = phi*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi + b + gamma*S(t))*D1(t),
    D2'(t) = phi*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    R'(t) = 0,
    k'(t) = 0,
    y1(t) = R(t),
    y2(t) = k(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kon kp y TT ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    C0'(t) = k(t) * (L1 - C0(t) - C1(t) - C2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao1) + phi(t)) * C0(t) + (b + gamma*S(t)) * C1(t),
    C1'(t) = phi(t)*C0(t) + (b + gamma*S(t))*C2(t) - ((1/tao1) + phi(t) + b + gamma*S(t))*C1(t),
    C2'(t) = phi(t)*C1(t) - ((1/tao1) + b + gamma*S(t))*C2(t),
    S'(t) = alpha*(C1(t)+D1(t))*(ST - S(t)) - beta*S(t),
    D0'(t) = k(t) * (L2 - D0(t) - D1(t) - D2(t)) * (R(t) - C0(t) - C1(t) - C2(t) - D0(t)- D1(t) - D2(t)) - ((1/tao2) + phi(t)) * D0(t) + (b + gamma*S(t)) * D1(t),
    D1'(t) = phi(t)*D0(t) + (b + gamma*S(t))*D2(t) - ((1/tao2) + phi(t) + b + gamma*S(t))*D1(t),
    D2'(t) = phi(t)*D1(t) - ((1/tao2) + b + gamma*S(t))*D2(t),
    phi'(t) = 0,
    k'(t) = 0,
    R'(t) = 0,
    y1(t) = R(t),
    y2(t) = phi(t),
    y3(t) = k(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))