# Dushek, O., & Van der Merwe, P. A. (2014). An induced rebinding model of antigen discrimination. 
#Trends in immunology, 35(4), 153-158.

using SIAN, Logging

# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff)*C1(t),
    y1(t) = C1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
