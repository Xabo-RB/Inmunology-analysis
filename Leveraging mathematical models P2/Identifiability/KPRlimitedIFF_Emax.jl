# Lever, M., Lim, H. S., Kruger, P., Nguyen, J., Trendel, N., Abu-Shah, E.,   & Dushek, O. (2016). 
# Architecture of a minimal signaling pathway explains the T-cell response to a 1 million-fold variation in antigen affinity and dose. 
# Proceedings of the National Academy of Sciences, 113(43), E6630-E6638.

using SIAN, Logging

# __________ todos los parámetros son distintos  ________________________________________________________
# -------------------  N = 2 -------------------
ode = @ODEmodel(
    #dPdt (pMHC) / dTdt (TCR) / dC0/dt (1º pMHC-TCR)
    L'(t) = -kon * L(t) * R(t) + koff * C0(t) + koff * C1(t) + koff * C2(t) + koff * C3(t) + koff * C4(t) + koff * C5(t) + koff * C6(t),
    R'(t) = - kon * L(t) * R(t) + koff*C0(t) + koff*C1(t) + koff*C2(t) + koff*C3(t) + koff*C4(t) + koff*C5(t) + koff*C6(t),
    C0'(t) = kon * L(t) * R(t) - (koff + kp)*C0(t),
    C1'(t) = kp*C0(t) - (koff + kp)*C1(t),
    C2'(t) = kp*C1(t) - (koff + kp)*C2(t),
    C3'(t) = kp*C2(t) - (koff + kp)*C3(t),
    C4'(t) = kp*C3(t) - (koff + kp)*C4(t),
    C5'(t) = kp*C4(t) - (koff + phi)*C5(t),
    C6'(t) = phi*C5(t) - koff*C6(t),
    Y'(t) = gyp * (YT - Y(t)) - gym * Y(t) + lambda * C5(t) * (YT - Y(t)),
    P'(t) = gpp * (PT - P(t)) - gpm * P(t) + Delta * Y(t) * (PT - P(t)) - mu * C5(t) * P(t), 
    y1(t) = (PT*(YT* Delta *(gyp + (((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) +
            (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu *
                 (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)*
                    mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta *
                     lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT*
                         Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT*
                         Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +
                       YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * mu *(-(gpm*gym*YT* Delta * lambda *
                  (YT* Delta *(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + 
                     ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )*
                 (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) *
                 (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* 
                      lambda ) + YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + 
                    (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) + 
                     (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))*
                  mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu *
              (-(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                  ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                 (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) -  
            sqrt((-4* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu *  
                    (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)*  
                       mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta *  
                        lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT*  
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *  
                   mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta *  
                       (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)*  
                          lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )*  
                    (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +   
                       ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu )))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*  
                (gpp + YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta *  
                    lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                 YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) +   
              ((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu *  
                     (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)*  
                        mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta *  
                         lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT*  
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT*  
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +   
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *  
                    mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta *  
                        (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)*  
                           lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )*  
                     (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                        ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                          ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                          ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp +  
                   YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                      lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +   
                   YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )))^  
               2)))/2) + gpp*(gym + gyp +  
         (((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2*  
                  mu *(-(gpm*(gpp + YT* Delta )* lambda ) +  
                  gpp*(gym + gyp)* mu  + gyp*YT* Delta * mu )*  
                 (gpm*gym*YT* Delta * lambda  - ((gym + gyp)*  
                      (gpp*(gym + gyp) + gyp*YT* Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*  
                      (gpp*(gym + gyp) + gyp*YT* Delta )* lambda  +  
                     ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp + YT* Delta )* lambda ^2)* mu )^ 
                  2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * mu *(-(gpm*gym*YT* Delta * lambda * 
                  (YT* Delta *(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp +  
                     ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                 (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) *  
                 (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* 
                      lambda ) + YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                    (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  
                     (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                  mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu *  
              (-(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                  ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                 (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) -  
            sqrt((-4* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu *  
                    (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)*  
                       mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta *  
                        lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT*  
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT*  
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *  
                   mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta *  
                       (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                          lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                    (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                       ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu )))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* 
                (gpp + YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                    lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                 YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) +  
              ((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                     (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                        mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                         lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                    mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                        (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                           lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                     (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                        ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                          ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                          ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp +  
                   YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                      lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                   YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )))^ 
               2)))/2)))/ 
     (YT* Delta *(gyp + (((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) +  
           (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu *  
                (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                   mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                    lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                        Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* Delta )* 
                      lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp + YT* Delta )* 
                      lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * mu * 
              (-(gpm*gym*YT* Delta * lambda *(YT* Delta *(gyp +  
                    ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ))) +  
               (gpp*(gym + gyp) + gyp*YT* Delta )*(gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon)  
                   +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) +  
                 YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                   (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon)  
                   +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                 mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu * 
             (-(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                 ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) -  
           sqrt((-4* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                   (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                      mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                       lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                         YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                  mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                         lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                   (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) *  
                   (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* 
                        lambda ) + YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  
                       (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                    mu )))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu *( 
                -(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                  (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) +  
             ((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                    (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                       mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                        lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                   mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                       (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                          lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                    (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                       ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp +  
                  YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                     lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                  YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )))^ 
              2)))/2) + (gym + gyp +  
        (((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda *((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                   mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                    lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                        Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* Delta )* 
                      lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp + YT* Delta )* 
                      lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * mu * 
              (-(gpm*gym*YT* Delta * lambda *(YT* Delta *(gyp +  
                    ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ))) +  
               (gpp*(gym + gyp) + gyp*YT* Delta )*(gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) +  
                 YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                   (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                 mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu * 
             (-(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                 ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) -  
           sqrt((-4* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                   (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                      mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                       lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                         YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                  mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                         lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                   (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon)  
                   +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* 
                        lambda ) + YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                    mu )))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu *( 
                -(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                  (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) +  
             ((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                    (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                       mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                        lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                   mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                       (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                          lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                    (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                       ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp +  
                  YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                     lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                  YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )))^ 
              2)))/2)*(gpm + gpp +  
        (((koff / (koff + phi)) * (kp/(koff+kp))^5)* mu *((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                   mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                    lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                        Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* Delta )* 
                      lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp + YT* Delta )* 
                      lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * mu * 
              (-(gpm*gym*YT* Delta * lambda *(YT* Delta *(gyp +  
                    ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ))) +  
               (gpp*(gym + gyp) + gyp*YT* Delta )*(gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) +  
                 YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                   (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                 mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu * 
             (-(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                 ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) -  
           sqrt((-4* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                   (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                      mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                       lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                         YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                  mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                         lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                   (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* 
                        lambda ) + YT* Delta *(gym*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) +  
                      (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* 
                    mu )))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp + YT* Delta )* lambda ^2* mu *( 
                -(gpm*gym*YT* Delta * lambda ) + gpp*(gym + gyp +  
                   ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  + YT* Delta *(gym*gyp +  
                  (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )) +  
             ((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)) + (sqrt(-(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*gym*YT* Delta * lambda ^2* mu * 
                    (-(gpm*(gpp + YT* Delta )* lambda ) + gpp*(gym + gyp)* 
                       mu  + gyp*YT* Delta * mu )*(gpm*gym*YT* Delta * 
                        lambda  - ((gym + gyp)*(gpp*(gym + gyp) + gyp*YT* 
                           Delta ) + 2*((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*(gpp*(gym + gyp) + gyp*YT* 
                           Delta )* lambda  + ((koff / (koff + phi)) * (kp/(koff+kp))^5)^2* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))*(gpp +  
                          YT* Delta )* lambda ^2)* mu )^2)) - ((koff / (koff + phi)) * (kp/(koff+kp))^5)* lambda * 
                   mu *(-(gpm*gym*YT* Delta * lambda *(YT* Delta * 
                       (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* 
                          lambda ))) + (gpp*(gym + gyp) + gyp*YT* Delta )* 
                    (gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gym + gyp +  
                       ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda ) + YT* Delta *(gym*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*(koff/kon)* lambda ) + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )*(gyp +  
                         ((koff / (koff + phi)) * (kp/(koff+kp))^5)*((koff/kon) +  (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t)))* lambda )))* mu ))/(((koff / (koff + phi)) * (kp/(koff+kp))^5)^2*(gpp +  
                  YT* Delta )* lambda ^2* mu *(-(gpm*gym*YT* Delta * 
                     lambda ) + gpp*(gym + gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2* mu  +  
                  YT* Delta *(gym*gyp + (gyp + ((koff / (koff + phi)) * (kp/(koff+kp))^5)* (R(t)+C0(t)+C1(t)+C2(t)+C3(t)+C4(t)+C5(t)+C6(t))* lambda )^2)* mu )))^ 
              2)))/2))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))
