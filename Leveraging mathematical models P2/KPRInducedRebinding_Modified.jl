# Dushek, O., & Van der Merwe, P. A. (2014). An induced rebinding model of antigen discrimination. 
#Trends in immunology, 35(4), 153-158.

using SIAN, Logging

# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    y1(t) = C2(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C2(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    y1(t) = C2(t),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO Kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t)* P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    y1(t) = C2(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kp'(t) = 0,
    y1(t) = C2(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ T_T ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    y1(t) = C2(t),
    y2(t) = C0(t) + C1(t) + C2(t) + T(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y = T(t)
# ===================================================================
# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------

ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    y1(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO Kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# ===================================================================
#               y = TT
# ===================================================================
# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    y1(t) = T(t) + C0(t) + C1(t) + C2(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO Kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon(t)* P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - kp*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp*C1(t) - (koff/(koff+rho2))*C2(t),
    kon'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + lambdaR*(koff/(koff+rho1))*C1(t) + lambdaR*(koff/(koff+rho2))*C2(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - kp(t)*C1(t) + (koff/(koff+rho1))*C1(t),
    C2'(t) = kp(t)*C1(t) - (koff/(koff+rho2))*C2(t),
    kp'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
