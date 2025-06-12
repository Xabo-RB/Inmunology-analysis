using StructuralIdentifiability

ode = StructuralIdentifiability.@ODEmodel(
    C_I'(t) = u(t) - etaa * (T_P(t) / (A + T_P(t))) * C_I(t) - mu_I * C_I(t),
    C_E'(t) = nuu * (T_P(t) / (A + T_P(t))) * C_I(t) + kappaa * (T_P(t) / (A + T_P(t))) * C_E(t) 
              - epsilonn * (1 - (T_P(t) / (A + T_P(t)))) * C_E(t) 
              + thetaa * (T_P(t) / (A + T_P(t))) * C_P(t) - mu_E * C_E(t),
    C_P'(t) = epsilonn * (1 - (T_P(t) / (A + T_P(t)))) * C_E(t) - thetaa * (T_P(t) / (A + T_P(t))) * C_P(t) - mu_P * C_P(t),
    T_P'(t) = rhoo * T_P(t) * (1 - (T_P(t) + T_N(t)) / K) - gammaa * (C_E(t) / (B + C_E(t))) * T_P(t),
    T_N'(t) = rhoo * T_N(t) * (1 - (T_P(t) + T_N(t)) / K) - g0 * gammaa * (C_E(t) / (B + C_E(t))) * T_N(t),
    M_i'(t) = sigma_M - (beta_B*(T_P(t)/(A + T_P(t)))*C_E(t) + beta_K*(C_E(t) / (B + C_E(t)))*(T_P(t) + g0*T_N(t))
                 + beta_C*(M_a(t) / (C + M_a(t)))*C_E(t))* M_i(t) - delta_M * M_i(t),
    M_a'(t) = (beta_B*(T_P(t)/(A + T_P(t)))*C_E(t) + beta_K*(C_E(t) / (B + C_E(t)))*(T_P(t) + g0*T_N(t))
                 + beta_C*(M_a(t) / (C + M_a(t)))*C_E)* M_i(t) - delta_M * M_a(t),
    IL_6'(t) = sigma_I + alphaa * M_a(t) - delta_I * IL_6(t),
    y1(t)= C_I(t),
    y2(t)= M_i(t) + M_a(t),
    y3(t)= IL_6(t)
)

println(assess_identifiability(ode))

using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*I(t)*S(t),
    E'(t) =  beta*I(t)*S(t) - alpha*E(t),
    I'(t) = alpha*E(t) - lambda*I(t),
    R'(t) = lambda*I(t),
    y1(t) = I(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10))