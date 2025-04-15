# Lever, M., Lim, H. S., Kruger, P., Nguyen, J., Trendel, N., Abu-Shah, E., ... & Dushek, O. (2016). 
# Architecture of a minimal signaling pathway explains the T-cell response to a 1 million-fold variation in antigen affinity and dose. 
# Proceedings of the National Academy of Sciences, 113(43), E6630-E6638.

using SIAN, Logging

# __________ todos los parámetros son distintos  ________________________________________________________
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
    y1(t) = (-1)*(-2 * lambda * nu - 2 * delta * lambda * nu - 3 * (gammaNeg/gammaPos) * lambda * nu - delta * (gammaNeg/gammaPos) * lambda * nu - (gammaNeg/gammaPos)^2 * lambda * nu - 
        2 * mu * nu - 3 * (gammaNeg/gammaPos) * mu * nu - (gammaNeg/gammaPos)^2 * mu * nu - lambda^2 * nu^2 - 2 * delta * lambda^2 * nu^2 - delta^2 * lambda^2 * nu^2 -
        (gammaNeg/gammaPos) * lambda^2 * nu^2 - delta * (gammaNeg/gammaPos) * lambda^2 * nu^2 - lambda * mu * nu^2 - delta * lambda * mu * nu^2 - (gammaNeg/gammaPos) * lambda * mu * nu^2 -
        delta * (gammaNeg/gammaPos) * lambda * mu * nu^2 - (nu * (2 + (gammaNeg/gammaPos) + lambda * nu + delta * lambda * nu) * 
        sqrt(lambda^2 + 2 * delta * lambda^2 + delta^2 * lambda^2 + 2 * (gammaNeg/gammaPos) * lambda^2 + 2 * delta * (gammaNeg/gammaPos) * lambda^2 + (gammaNeg/gammaPos)^2 * lambda^2 - 2 * 
        lambda * mu - 2 * delta * lambda * mu - 4 * (gammaNeg/gammaPos) * lambda * mu + 2 * delta * (gammaNeg/gammaPos) * lambda * mu - 2 * (gammaNeg/gammaPos)^2 * lambda * mu + mu^2 + 2 * 
        (gammaNeg/gammaPos) * mu^2 + (gammaNeg/gammaPos)^2 * mu^2)))/(2 * (2 * lambda * mu * nu^2 + (gammaNeg/gammaPos) * lambda * mu * nu^2 + lambda^2 * mu * nu^3 + delta * lambda^2 * mu * nu^3))
        * (R(t) + C0(t) + C1(t) + C2(t)) + sqrt(-(koff/kon)^2 - 2 * (-2 * lambda * nu - 2 * delta * lambda * nu - 3 * (gammaNeg/gammaPos) * lambda * nu - delta * (gammaNeg/gammaPos) * lambda * nu - (gammaNeg/gammaPos)^2 * lambda * nu - 
        2 * mu * nu - 3 * (gammaNeg/gammaPos) * mu * nu - (gammaNeg/gammaPos)^2 * mu * nu - lambda^2 * nu^2 - 2 * delta * lambda^2 * nu^2 - delta^2 * lambda^2 * nu^2 -
        (gammaNeg/gammaPos) * lambda^2 * nu^2 - delta * (gammaNeg/gammaPos) * lambda^2 * nu^2 - lambda * mu * nu^2 - delta * lambda * mu * nu^2 - (gammaNeg/gammaPos) * lambda * mu * nu^2 -
        delta * (gammaNeg/gammaPos) * lambda * mu * nu^2 - (nu * (2 + (gammaNeg/gammaPos) + lambda * nu + delta * lambda * nu) * 
        sqrt(lambda^2 + 2 * delta * lambda^2 + delta^2 * lambda^2 + 2 * (gammaNeg/gammaPos) * lambda^2 + 2 * delta * (gammaNeg/gammaPos) * lambda^2 + (gammaNeg/gammaPos)^2 * lambda^2 - 2 * 
        lambda * mu - 2 * delta * lambda * mu - 4 * (gammaNeg/gammaPos) * lambda * mu + 2 * delta * (gammaNeg/gammaPos) * lambda * mu - 2 * (gammaNeg/gammaPos)^2 * lambda * mu + mu^2 + 2 * 
        (gammaNeg/gammaPos) * mu^2 + (gammaNeg/gammaPos)^2 * mu^2)))/(2 * (2 * lambda * mu * nu^2 + (gammaNeg/gammaPos) * lambda * mu * nu^2 + lambda^2 * mu * nu^3 + delta * lambda^2 * mu * nu^3))
        * (koff/kon) * (R(t) + C0(t) + C1(t) + C2(t)) - 2 * (-2 * lambda * nu - 2 * delta * lambda * nu - 3 * (gammaNeg/gammaPos) * lambda * nu - delta * (gammaNeg/gammaPos) * lambda * nu - (gammaNeg/gammaPos)^2 * lambda * nu - 
        2 * mu * nu - 3 * (gammaNeg/gammaPos) * mu * nu - (gammaNeg/gammaPos)^2 * mu * nu - lambda^2 * nu^2 - 2 * delta * lambda^2 * nu^2 - delta^2 * lambda^2 * nu^2 -
        (gammaNeg/gammaPos) * lambda^2 * nu^2 - delta * (gammaNeg/gammaPos) * lambda^2 * nu^2 - lambda * mu * nu^2 - delta * lambda * mu * nu^2 - (gammaNeg/gammaPos) * lambda * mu * nu^2 -
        delta * (gammaNeg/gammaPos) * lambda * mu * nu^2 - (nu * (2 + (gammaNeg/gammaPos) + lambda * nu + delta * lambda * nu) * 
        sqrt(lambda^2 + 2 * delta * lambda^2 + delta^2 * lambda^2 + 2 * (gammaNeg/gammaPos) * lambda^2 + 2 * delta * (gammaNeg/gammaPos) * lambda^2 + (gammaNeg/gammaPos)^2 * lambda^2 - 2 * 
        lambda * mu - 2 * delta * lambda * mu - 4 * (gammaNeg/gammaPos) * lambda * mu + 2 * delta * (gammaNeg/gammaPos) * lambda * mu - 2 * (gammaNeg/gammaPos)^2 * lambda * mu + mu^2 + 2 * 
        (gammaNeg/gammaPos) * mu^2 + (gammaNeg/gammaPos)^2 * mu^2)))/(2 * (2 * lambda * mu * nu^2 + (gammaNeg/gammaPos) * lambda * mu * nu^2 + lambda^2 * mu * nu^3 + delta * lambda^2 * mu * nu^3))
        * (R(t) + C0(t) + C1(t) + C2(t)) + 5 * ((-2 * lambda * nu - 2 * delta * lambda * nu - 3 * (gammaNeg/gammaPos) * lambda * nu - delta * (gammaNeg/gammaPos) * lambda * nu - (gammaNeg/gammaPos)^2 * lambda * nu - 
        2 * mu * nu - 3 * (gammaNeg/gammaPos) * mu * nu - (gammaNeg/gammaPos)^2 * mu * nu - lambda^2 * nu^2 - 2 * delta * lambda^2 * nu^2 - delta^2 * lambda^2 * nu^2 -
        (gammaNeg/gammaPos) * lambda^2 * nu^2 - delta * (gammaNeg/gammaPos) * lambda^2 * nu^2 - lambda * mu * nu^2 - delta * lambda * mu * nu^2 - (gammaNeg/gammaPos) * lambda * mu * nu^2 -
        delta * (gammaNeg/gammaPos) * lambda * mu * nu^2 - (nu * (2 + (gammaNeg/gammaPos) + lambda * nu + delta * lambda * nu) * 
        sqrt(lambda^2 + 2 * delta * lambda^2 + delta^2 * lambda^2 + 2 * (gammaNeg/gammaPos) * lambda^2 + 2 * delta * (gammaNeg/gammaPos) * lambda^2 + (gammaNeg/gammaPos)^2 * lambda^2 - 2 * 
        lambda * mu - 2 * delta * lambda * mu - 4 * (gammaNeg/gammaPos) * lambda * mu + 2 * delta * (gammaNeg/gammaPos) * lambda * mu - 2 * (gammaNeg/gammaPos)^2 * lambda * mu + mu^2 + 2 * 
        (gammaNeg/gammaPos) * mu^2 + (gammaNeg/gammaPos)^2 * mu^2)))/(2 * (2 * lambda * mu * nu^2 + (gammaNeg/gammaPos) * lambda * mu * nu^2 + lambda^2 * mu * nu^3 + delta * lambda^2 * mu * nu^3)))^2
        * (R(t) + C0(t) + C1(t) + C2(t))^2)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ GtHETA  ________________________________________________________
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
    y1(t) = -gtheta * (R(t) + C0(t) + C1(t) + C2(t)) + (-(koff/kon)^2 - 2 * g_theta * (koff/kon) * (R(t) + C0(t) + C1(t) + C2(t)) - 2 * gtheta * (R(t) + C0(t) + C1(t) + C2(t)) + 5 * gtheta^2 * (R(t) + C0(t) + C1(t) + C2(t))^2)^(1/2)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

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
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = R(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ T(t)  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * L(t) * R(t)  - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi*kp)*C1(t),
    C2'(t) = phi*kp*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = R(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

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
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

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
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ TT(t) ________________________________________________________
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
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = R(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Conociendo todos  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * L(t) * R(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi*kp(t))*C1(t),
    C2'(t) = phi*kp(t)*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = PT * (2 + (gammaNeg/gammaPos) + (lambda + delta * lambda) * mu) / (lambda * mu * nu^2 + (lambda + lambda * (gammaNeg/gammaPos) + mu + (gammaNeg/gammaPos) * mu + delta * lambda) * nu + (gammaNeg/gammaPos)^2 + 2 * (gammaNeg/gammaPos) + delta + 1),
    y2(t) = kon(t),
    y3(t) = kp(t),
    y4(t) = R(t) + C0(t) + C1(t) + C2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon(t) * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon(t) * L(t) * R(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + phi*kp)*C1(t),
    C2'(t) = phi*kp*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kon'(t) = 0,
    y1(t) = R(t),
    y2(t) = kon(t)
)

# __________ kp  ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t),
    C0'(t) = kon * L(t) * R(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + phi*kp(t))*C1(t),
    C2'(t) = phi*kp(t)*C1(t) - koff*C2(t),
    Y'(t) = gammaPos * (YT - Y(t)) - gammaNeg*Y(t) + lambda*C1(t)*(YT - Y(t)),
    P'(t) = gammaPos * (PT - P(t)) - gammaNeg*P(t) + delta*Y(t)*(PT - P(t)) - mu*C1(t)*P(t), 
    kp'(t) = 0,
    y1(t) = R(t),
    y2(t) = kp(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
