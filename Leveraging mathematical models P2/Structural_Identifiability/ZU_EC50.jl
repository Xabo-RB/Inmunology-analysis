#Goldbeter, A., & Koshland Jr, D. E. (1981). An amplified sensitivity arising from covalent modification in biological systems. 
# Proceedings of the National Academy of Sciences, 78(11), 6840-6844.

using SIAN, Logging

# __________ Nada _______________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1 * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1 * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    y1(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1*(k3 + kmenos2))) - 3*k1*w*(k3 + kmenos2)^2)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ T(t) ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1 * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1 * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    y1(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1*(k3 + kmenos2))) - 3*k1*w*(k3 + kmenos2)^2),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ kon ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1(t) * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1(t) * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    k1'(t) = 0,
    y1(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1(t)*(k3 + kmenos2))) - 3*k1(t)*w*(k3 + kmenos2)^2),
    y2(t) = k1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Tt ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1 * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1 * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    y1(t) = T(t) + Tp(t) + C0(t) + D(t),
    y2(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1*(k3 + kmenos2))) - 3*k1*w*(k3 + kmenos2)^2)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))


# __________ T y kon ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1(t) * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1(t) * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    k1'(t) = 0,
    y1(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1(t)*(k3 + kmenos2))) - 3*k1(t)*w*(k3 + kmenos2)^2),
    y2(t) = k1(t),
    y3(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ Tt y kon ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1(t) * T(t) * L(t) + k3 * D(t) + kmenos1 * C0(t),
    L'(t) = -k1(t) * T(t) * L(t) + w * C0(t) + kmenos1 * C0(t),
    C0'(t) = k1(t) * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    k1'(t) = 0,
    y1(t) = (2*k2*k3*(Q(t)+D(t))*(k3 + kmenos2 - k2*(T(t) + Tp(t) + C0(t) + D(t)))*(kmenos1+w))/ ( k2*(T(t) + Tp(t) + C0(t) + D(t))*(2*kmenos1 + w*(2 + k1(t)*(k3 + kmenos2))) - 3*k1(t)*w*(k3 + kmenos2)^2),
    y2(t) = k1(t),
    y3(t) = T(t) + Tp(t) + C0(t) + D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))