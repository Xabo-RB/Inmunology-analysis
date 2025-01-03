# Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of 
# ligand discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.


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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2 # T_T \cdot \left(\dfrac{k_{-2} + k_3}{k_2}\right)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T/(t) ________________________________________________________
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T_T ________________________________________________________
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
    y2(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ KON ________________________________________________________
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
    y2(t) = w(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp T_T ________________________________________________________
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)

)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp T ________________________________________________________
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
    y1(t) = (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t))*kmenos2*k3/k2,
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = T(t)

)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# ===================================================================
#               y =  EC50
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w)^2 - w * (kmenos1 + 2*w) * (w / (kmenos1 + w))^3)) / (k1 * w * (w / (kmenos1 + w))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1 * (Q(t) + D(t)) * (kmenos1 + w))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T/(t) ________________________________________________________
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w)^2 - w * (kmenos1 + 2*w) * (w / (kmenos1 + w))^3)) / (k1 * w * (w / (kmenos1 + w))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1 * (Q(t) + D(t)) * (kmenos1 + w)),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T_T ________________________________________________________
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w)^2 - w * (kmenos1 + 2*w) * (w / (kmenos1 + w))^3)) / (k1 * w * (w / (kmenos1 + w))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1 * (Q(t) + D(t)) * (kmenos1 + w)),
    y2(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ KON ________________________________________________________
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w)^2 - w * (kmenos1 + 2*w) * (w / (kmenos1 + w))^3)) / (k1(t) * w * (w / (kmenos1 + w))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1(t) * (Q(t) + D(t)) * (kmenos1 + w)),
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w(t))^2 - w(t) * (kmenos1 + 2*w(t)) * (w(t) / (kmenos1 + w(t)))^3)) / (k1 * w(t) * (w(t) / (kmenos1 + w(t)))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1 * (Q(t) + D(t)) * (kmenos1 + w)),
    y2(t) = w(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp T_T ________________________________________________________
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w(t))^2 - w(t) * (kmenos1 + 2*w(t)) * (w(t) / (kmenos1 + w(t)))^3)) / (k1(t) * w(t) * (w(t) / (kmenos1 + w(t)))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1(t) * (Q(t) + D(t)) * (kmenos1 + w)),
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon kp T ________________________________________________________
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
    y1(t) = (2 * k2 * k3 * (Q(t) + D(t)) * ((kmenos1 + w(t))^2 - w(t) * (kmenos1 + 2*w(t)) * (w(t) / (kmenos1 + w(t)))^3)) / (k1(t) * w(t) * (w(t) / (kmenos1 + w(t)))^3 * (-2 * k2 * k3 * (Q(t) + D(t)) + 2 * k2 * kmenos1 * ((T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t)) + 2) + kmenos1 * (T(t) + C0(t) + C1(t) + C2(t) + Tp(t) + D(t) + 2) * (k3 + kmenos2)) + 2 * k2 * k3 * k1(t) * (Q(t) + D(t)) * (kmenos1 + w)),
    y2(t) = w(t),
    y3(t) = k1(t),
    y4(t) = T(t)

)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

