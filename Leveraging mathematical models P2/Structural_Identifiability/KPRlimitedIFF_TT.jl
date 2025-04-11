# Lever, M., Lim, H. S., Kruger, P., Nguyen, J., Trendel, N., Abu-Shah, E., ... & Dushek, O. (2016). 
# Architecture of a minimal signaling pathway explains the T-cell response to a 1 million-fold variation in antigen affinity and dose. 
# Proceedings of the National Academy of Sciences, 113(43), E6630-E6638.

using SIAN, Logging

# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * L(t) * R(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi*kp)*C1(t),
    C2'(t) = phi*kp*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    y1(t) = R(t) + C0(t) + C1(t) + C2(t),
)

# __________ Conociendo todos  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * L(t) * R(t)  - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi*kp(t))*C1(t),
    C2'(t) = phi*kp(t)*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = P(t),
    y2(t) = R(t) + C0(t) + C1(t) + C2(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# __________ Conociendo todos y no respuesta ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * L(t) * R(t)  - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi*kp(t))*C1(t),
    C2'(t) = phi*kp(t)*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y2(t) = R(t) + C0(t) + C1(t) + C2(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# __________ Respuesta(t)  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * L(t) * R(t)  - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi*kp)*C1(t),
    C2'(t) = phi*kp*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    y1(t) = R(t) + C0(t) + C1(t) + C2(t),
    y2(t) = P(t)
)

# __________ kon  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * L(t) * R(t)  - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi*kp)*C1(t),
    C2'(t) = phi*kp*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kon'(t) = 0,
    y1(t) = R(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kon(t)
)

# __________ kp  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * L(t) * R(t)  - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi*kp(t))*C1(t),
    C2'(t) = phi*kp(t)*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kp'(t) = 0,
    y1(t) = R(t) + C0(t) + C1(t) + C2(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

