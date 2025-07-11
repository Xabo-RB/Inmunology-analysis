# Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of ligand 
# discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.

using SIAN, Logging

# __________ nada ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1 * T(t) * X(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * X(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = k1 * (kmenos1 + k3) * w^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w)^3 + k1 * (kmenos1 + k3) * w^3)
)
    #y1(t) = ((C0(t)+C1(t)+C2(t)) * (k2 + k3) * kmenos1 * w * (w / (kmenos1 + w))^2)/(k2 * (k2 + k3 - kmenos2) * (kmenos1 + w - w * (w / (kmenos1 + w))^2) * X(t))


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# -------------------  N = Tt -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1 * T(t) * X(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * X(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) =  k1 * (kmenos1 + k3) * w^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w)^3 + k1 * (kmenos1 + k3) * w^3),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# -------------------  kon -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1(t) * T(t) * X(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * X(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    k1'(t) = 0,
    y1(t) = k1(t) * (kmenos1 + k3) * w^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w)^3 + k1(t) * (kmenos1 + k3) * w^3),
    y2(t) = k1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1 * T(t) * X(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * X(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    w'(t) = 0,
    y1(t) = w(t),
    y2(t) = k1 * (kmenos1 + k3) * w(t)^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w(t))^3 + k1 * (kmenos1 + k3) * w(t)^3)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ todos con TT ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1(t) * T(t) * X(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * X(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    k1'(t) = 0,
    w'(t) = 0,
    y1(t) = w(t),
    y1(t) = k1(t)*(kmenos1 + k3) * w(t)^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w(t))^3 + k1(t)* (kmenos1 + k3) * w(t)^3),
    y3(t) = k1(t),    
    y4(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t),

)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ todos con T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1(t) * T(t) * X(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * X(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    k1'(t) = 0,
    w'(t) = 0,
    y1(t) = w(t),
    y2(t) = k1(t)*(kmenos1 + k3) * w(t)^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w(t))^3 + k1(t)* (kmenos1 + k3) * w(t)^3),
    y3(t) = k1(t),    
    y4(t) = T(t)

)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# -------------------  N = T -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * X(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    X'(t) = -k1 * T(t) * X(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t) + k3*D(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*X(t) + kmenos2*D(t),
    D'(t) = k2 * Tp(t) * X(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * X(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = k1 * (kmenos1 + k3) * w^3*(T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) /(k2 * k3 * (kmenos1 + w)^3 + k1 * (kmenos1 + k3) * w^3),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
