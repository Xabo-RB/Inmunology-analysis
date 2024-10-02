#Gálvez, J., Galvez, J. J., & Garcia-Penarrubia, P. (2016). 
#TCR/pMHC interaction: phenotypic model for an unsolved enigma. 
#Frontiers in Immunology, 7, 228066.

using SIAN, Logging


# __________ nada ________________________________________________________
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t) + C4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ R(t) ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t) + (4/(1 + 3*r))*koff*C3(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - ((3/(1 + 2*r))*koff + (rp^2)*kp)*C2(t),
    C3'(t) = (3/(1 + 2*r))*kp*C2(t) - (4/(1 + 3*r))*koff*C3(t),
    y1(t) = C3(t),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t) + C3(t)
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kp(t)
)

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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t)
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
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y =  Emax 
# ===================================================================

# __________ nada ________________________________________________________
# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t)
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
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kp(t)
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
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ TT ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS TT ________________________________________________________
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
    y4(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS T ________________________________________________________
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
    y4(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y =  EC50 
# ===================================================================

# __________ nada ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = (kon*(24*koff^3 + kp^3*(1+6*r+11*r^2+6*r^3)*rp^3 + 4*koff*kp^2*(1+3*r+2*r^2)*rp*(2+rp) + 4*koff^2*kp*(3+5*rp+r*(3+7*rp)))) / (4*koff*(koff+kp)*(2*koff+kp*(1+r)*rp)*(3*koff+kp*(1+2*r)*rp)) + (T(t) + C0(t) + C1(t) + C2(t))/2
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t)
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
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kp(t)
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
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = kon(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ TT ________________________________________________________
# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + (2/(1+r))*koff*C1(t) + (3/(1+2*r))*koff*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - ((2/(1+r))*koff + rp*kp)*C1(t),
    C2'(t) = rp*kp*C1(t) - (3/(1 + 2*r))*koff*C2(t),
    y1(t) = ((1 + (4 * koff * (6 * koff^2 + kp^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t)),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS TT ________________________________________________________
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
    y4(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________  TODOS T ________________________________________________________
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
    y4(t) = ((1 + (4 * koff * (6 * koff^2 + kp(t)^2 * (1 + 3 * r + 2 * r^2) * rp^3 * (2 + rp) + koff * kp(t) * rp * (3 + 5 * rp + r * (3 + 7 * rp)))) / (kp(t)^3 * (1 + r) * (1 + 2 * r) * (1 + 3 * r) * rp^5)) ) * (T(t) + C0(t) + C1(t) + C2(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))