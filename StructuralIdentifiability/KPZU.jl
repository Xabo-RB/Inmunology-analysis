# Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of ligand 
# discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.

using SIAN, Logging


# __________ nada ________________________________________________________
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    T'(t) = -k1 * T(t) * P(t) + k3 * D(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    P'(t) = -k1 * T(t) * P(t) + w*C2(t) + kmenos1 * C0(t) + kmenos1*C1(t) + kmenos1*C2(t),
    
    C0'(t) = k1 * T(t) * L(t) - (kmenos1 + w) * C0(t),
    D'(t) = k2 * Tp(t) * Q(t) - (kmenos2 + k3) * D(t),
    Tp'(t) = -k2 * Tp(t) * Q(t) +kmenos2 * D(t) + w * C0(t),
    Q'(t) = -k2 * Tp(t) * Q(t) + kmenos2 * D(t) + k3 * D(t),
    y1(t) = Tp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))