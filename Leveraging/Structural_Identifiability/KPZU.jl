# Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of ligand 
# discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.

using SIAN, Logging


# __________ nada ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = Tp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Tt ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = Tp(t),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    k1'(t) = 0,
    y1(t) = Tp(t),
    y2(t) = k1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    w'(t) = 0,
    y1(t) = Tp(t),
    y2(t) = w(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ todos ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    k1'(t) = 0,
    w'(t) = 0,
    y1(t) = Tp(t),
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)

)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y =  T(t)
# ===================================================================


# __________ nada ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ R ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = T(t),
    y2(t) = Tp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    k1'(t) = 0,
    y1(t) = T(t),
    y2(t) = k1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    w'(t) = 0,
    y1(t) = T(t),
    y2(t) = w(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ todos ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    k1'(t) = 0,
    w'(t) = 0,
    y1(t) = T(t),
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = R(t)

)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y =  Tt(t)
# ===================================================================


# __________ nada ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w*C1(t) + w*C0(t),
    C2'(t) = - kmenos1*C2(t) - w*C2(t) + w*C1(t),
    k1'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t),
    y2(t) = k1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1 * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    w'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t),
    y2(t) = w(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ todos ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1(t) * T(t) * P(t) + w(t)*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    Tp'(t) = w(t)*C2(t) - k2*Tp(t)*Q(t) + kmenos2*D(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    C0'(t) = k1(t) * T(t) * P(t) - (kmenos1 + w(t)) * C0(t),
    C1'(t) = - kmenos1*C1(t) - w(t)*C1(t) + w(t)*C0(t),
    C2'(t) = - kmenos1*C2(t) - w(t)*C2(t) + w(t)*C1(t),
    k1'(t) = 0,
    w'(t) = 0,
    y1(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t),
    y2(t) = w(t),
    y3(t) = k1(t)

)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))