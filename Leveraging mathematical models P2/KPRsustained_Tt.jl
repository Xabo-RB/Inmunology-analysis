#Coombs, D., Kalergis, A. M., Nathenson, S. G., Wofsy, C., & Goldstein, B. (2002). 
#Activated TCRs remain marked for internalization after dissociation from pMHC. 
#Nature immunology_, _3_(10), 926-931.

#González, P. A., Carreño, L. J., Coombs, D., Mora, J. E., Palmieri, E., Goldstein, B., ... & Kalergis, A. M. (2005). 
#T cell receptor binding kinetics required for T cell activation depend on the density of cognate ligand on the antigen-presenting cell. 
#_Proceedings of the National Academy of Sciences_, _102_(13), 4824-4829.

using SIAN, Logging

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) - kon * P(t) * Tast(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + lambda * Tast(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - koff*C3(t) + kon*P(t)*Tast(t),
    Tast'(t) = koff*C3(t) - kon*P(t)*Tast(t) - lambda*Tast(t), 
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + Tast(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo TODOS ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) - kon(t) * P(t) * Tast(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + lambda * Tast(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - koff*C3(t) + kon(t)*P(t)*Tast(t),
    Tast'(t) = koff*C3(t) - kon(t)*P(t)*Tast(t) - lambda*Tast(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + Tast(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = C3(t) + Tast(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo R(t) ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) - kon * P(t) * Tast(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + lambda * Tast(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - koff*C3(t) + kon*P(t)*Tast(t),
    Tast'(t) = koff*C3(t) - kon*P(t)*Tast(t) - lambda*Tast(t), 
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + Tast(t),
    y2(t) = C3(t) + Tast(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo KON ________________________________________________________
# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) - kon(t) * P(t) * Tast(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + lambda * Tast(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - koff*C4(t) + kon(t)*P(t)*Tast(t),
    Tast'(t) = koff*C4(t) - kon(t)*P(t)*Tast(t) - lambda*Tast(t), 
    kon'(t) = 0, 
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + C4(t) +  Tast(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo kp ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) - kon * P(t) * Tast(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + lambda * Tast(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t),
    C4'(t) = kp(t)*C3(t) - koff*C4(t) + kon*P(t)*Tast(t),
    Tast'(t) = koff*C4(t) - kon*P(t)*Tast(t) - lambda*Tast(t), 
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + C4(t) +  Tast(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ conociendo parametros ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) - kon(t) * P(t) * Tast(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + lambda * Tast(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t),
    C3'(t) = kp(t)*C2(t) - koff*C3(t) + kon(t)*P(t)*Tast(t),
    Tast'(t) = koff*C3(t) - kon(t)*P(t)*Tast(t) - lambda*Tast(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + Tast(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
