# McKeithan, T. W. (1995). Kinetic proofreading in T-cell receptor signal transduction. 
#_Proceedings of the national academy of sciences_, _92_(11), 5042-5046.

using SIAN, Logging

# -------------------  N = 1 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    T'(t) = - kon * P(t) * T(t) + koff*C0(t) + koff*C1(t),
    C0'(t) = - kon * P * T - koff*C0(t),
    y1(t) = C0(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))