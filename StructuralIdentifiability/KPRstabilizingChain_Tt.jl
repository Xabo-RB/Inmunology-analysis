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
    y1(t) = C4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

