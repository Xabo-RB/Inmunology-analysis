# Dushek, O., & Van der Merwe, P. A. (2014). An induced rebinding model of antigen discrimination. 
#Trends in immunology, 35(4), 153-158.

using SIAN, Logging

# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    y1(t) = C2a(t) + C2(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    y1(t) = C3a(t) + C3(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    y1(t) = C4a(t) + C4(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t) + rho4*C4a(t),
    C5'(t) = kp*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp)*C4a(t), 
    C5a'(t) = kp*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    y1(t) = C5a(t) + C5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C2a(t) + C2(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C3a(t) + C3(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C4a(t) + C4(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t) + rho4*C4a(t),
    C5'(t) = kp(t)*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp(t))*C4a(t), 
    C5a'(t) = kp(t)*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = C5a(t) + C5(t),
    y2(t) = T(t),
    y3(t) = kon(t),
    y4(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO T(t) ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    y1(t) = C2a(t) + C2(t),
    y2(t) = T(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    y1(t) = C3a(t) + C3(t),
    y2(t) = T(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    y1(t) = C4a(t) + C4(t),
    y2(t) = T(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t) + rho4*C4a(t),
    C5'(t) = kp*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp)*C4a(t), 
    C5a'(t) = kp*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    y1(t) = C5a(t) + C5(t),
    y2(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO Kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kon'(t) = 0,
    y1(t) = C2a(t) + C2(t),
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kon'(t) = 0,
    y1(t) = C3a(t) + C3(t),
    y2(t) = kon(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kon'(t) = 0,
    y1(t) = C4a(t) + C4(t),
    y2(t) = kon(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t) + rho4*C4a(t),
    C5'(t) = kp*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp)*C4a(t), 
    C5a'(t) = kp*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kon'(t) = 0,
    y1(t) = C5a(t) + C5(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon* P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kp'(t) = 0,
    y1(t) = C2a(t) + C2(t),
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kp'(t) = 0,
    y1(t) = C3a(t) + C3(t),
    y2(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kp'(t) = 0,
    y1(t) = C4a(t) + C4(t),
    y2(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t) + rho4*C4a(t),
    C5'(t) = kp(t)*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp(t))*C4a(t), 
    C5a'(t) = kp(t)*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kp'(t) = 0,
    y1(t) = C5a(t) + C5(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))



# ===================================================================
#               y = T(t)
# ===================================================================
# __________ SIN CONOCER NINGÚN PARÁMETRO ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t) + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t) + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    y1(t) = T(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    y1(t) = T(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    y1(t) = T(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t) + rho4*C4a(t),
    C5'(t) = kp*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp)*C4a(t), 
    C5a'(t) = kp*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    y1(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO TODOS ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t) + rho4*C4a(t),
    C5'(t) = kp(t)*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp(t))*C4a(t), 
    C5a'(t) = kp(t)*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kon'(t) = 0,
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t),
    y3(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO Kon ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon(t)*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon(t) * P(t) * T(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t) + rho1*C1a(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t) + rho2*C2a(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t) + rho3*C3a(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t) + rho4*C4a(t),
    C5'(t) = kp*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp)*C1a(t),
    C2a'(t) = kp*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp)*C2a(t),    
    C3a'(t) = kp*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp)*C3a(t), 
    C4a'(t) = kp*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp)*C4a(t), 
    C5a'(t) = kp*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kon'(t) = 0,
    y1(t) = T(t),
    y2(t) = kon(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

# __________ CONOCIENDO kp ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t),
    C0'(t) = kon* P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff)*C2(t) + rho2*C2a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR)*C2a(t),  
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 3 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff)*C3(t) + rho3*C3a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR)*C3a(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 4 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff)*C4(t) + rho4*C4a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR)*C4a(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

# -------------------  N = 5 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    P'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    T'(t) = - kon*P(t)*T(t) + koff*C0(t)  + koff*C1a(t) + koff*C2a(t) + koff*C3a(t) + koff*C4a(t) + koff*C5a(t),
    C0'(t) = kon * P(t) * T(t) - (koff + kp(t))*C0(t),
    C1'(t) = kp(t)*C0(t) - (koff + kp(t))*C1(t) + rho1*C1a(t),
    C2'(t) = kp(t)*C1(t) - (koff + kp(t))*C2(t) + rho2*C2a(t),
    C3'(t) = kp(t)*C2(t) - (koff + kp(t))*C3(t) + rho3*C3a(t),
    C4'(t) = kp(t)*C3(t) - (koff + kp(t))*C4(t) + rho4*C4a(t),
    C5'(t) = kp(t)*C4(t) - (koff)*C5(t) + rho5*C5a(t),
    C1a'(t) = koff*C1(t) - (rho1 + lambdaR + kp(t))*C1a(t),
    C2a'(t) = kp(t)*C1a(t) + koff*C2(t) - (rho2 + lambdaR + kp(t))*C2a(t),    
    C3a'(t) = kp(t)*C2a(t) + koff*C3(t) - (rho3 + lambdaR + kp(t))*C3a(t), 
    C4a'(t) = kp(t)*C3a(t) + koff*C4(t) - (rho4 + lambdaR + kp(t))*C4a(t), 
    C5a'(t) = kp(t)*C4a(t) + koff*C5(t) - (rho5 + lambdaR)*C5a(t),
    kp'(t) = 0,
    y1(t) = T(t),
    y2(t) = kp(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
