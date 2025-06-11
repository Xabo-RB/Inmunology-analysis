#Gálvez, J., Galvez, J. J., & Garcia-Penarrubia, P. (2016). 
#TCR/pMHC interaction: phenotypic model for an unsolved enigma. 
#Frontiers in Immunology, 7, 228066.

using SIAN, Logging


# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp*rp*(2*r+1)) * ( 1 + (2*koff)/(kp*(r+1)) + rp ))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp*rp*(2*r+1)) * ( 1 + (2*koff)/(kp*(r+1)) + rp )),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + rp*kp(t))*C1(t),
    C2'(t) = rp*kp(t)*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp(t)*rp*(2*r+1)) * ( 1 + (2*koff)/(kp(t)*(r+1)) + rp ))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kon'(t) = 0,
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp*rp*(2*r+1)) * ( 1 + (2*koff)/(kp*(r+1)) + rp )),
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
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + rp*kp(t))*C1(t),
    C2'(t) = rp*kp(t)*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kp'(t) = 0,
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp(t)*rp*(2*r+1)) * ( 1 + (2*koff)/(kp(t)*(r+1)) + rp )),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T_T ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp*rp*(2*r+1)) * ( 1 + (2*koff)/(kp*(r+1)) + rp )),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS CON TT________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon(t) * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - ((2/(1+r))*koff + rp*kp(t))*C1(t),
    C2'(t) = rp*kp(t)*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = (T(t) + C0(t) + C1(t) + C2(t)) * (1 + (3*koff)/(kp(t)*rp*(2*r+1)) * ( 1 + (2*koff)/(kp(t)*(r+1)) + rp ))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))