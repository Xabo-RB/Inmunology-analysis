
import sympy
from sympy import symbols, sqrt, simplify, solve, diff
from sympy import factor, radsimp, cancel, latex

# 1) Definir variables simbólicas
koff, kon, kp = symbols("koff kon kp", positive=True)
LT, TT, PT = symbols("LT TT PT", positive=True)
alpha = symbols("alpha", positive=True)
phi = symbols("phi", positive=True)   # Phi -> phi
delta = symbols("delta", positive=True)  # Delta -> delta
lmbd = symbols("lambda", positive=True)  # Lambda -> lmbd (evito 'lambda' que es palabra reservada)
mu = symbols("mu", positive=True)   
Gamma_poff, Gamma_pon = symbols("Gamma_poff Gamma_pon", positive=True)
Gamma_yoff, Gamma_yon = symbols("Gamma_yoff Gamma_yon", positive=True)
YT = symbols("YT", positive=True)

# 2) Expresión tal como la pegaste, con ligeros ajustes de sintaxis.
#    Notar que en Mathematica se usan backslash y a veces ^ para exponentes, 
#    mientras que en Python/Sympy es **. También cambio \[Gamma]yoff a Gamma_yoff, etc.
expr = (4*PT*(koff + phi)*(koff*Gamma_pon*(Gamma_yoff + Gamma_yon)
          + koff*YT*Gamma_yon*delta
          - ( kp*(koff + kon*(LT + TT) - sqrt(-4*kon**2*LT*TT + (koff + kon*(LT + TT))**2))
               *(-1 + alpha)*alpha**4*(Gamma_pon + YT*delta)*lmbd )/(2*kon)
          + Gamma_pon*(Gamma_yoff + Gamma_yon)*phi
          + YT*Gamma_yon*delta*phi ))
expr /= ( 4*koff**2*((Gamma_poff + Gamma_pon)*(Gamma_yoff + Gamma_yon) + YT*Gamma_yon*delta)
          + ( kp**2*(koff + kon*(LT + TT)
                     - sqrt(-4*kon**2*LT*TT + (koff + kon*(LT + TT))**2))**2
               *(-1 + alpha)**2*alpha**8*lmbd*mu )/kon**2
          - (1/kon)*2*koff*kp*(koff + kon*(LT + TT)
             - sqrt(-4*kon**2*LT*TT + (koff + kon*(LT + TT))**2))
             *(-1 + alpha)*alpha**4
             *((Gamma_poff + Gamma_pon + YT*delta)*lmbd + (Gamma_yoff + Gamma_yon)*mu)
          + 8*koff*((Gamma_poff + Gamma_pon)*(Gamma_yoff + Gamma_yon) + YT*Gamma_yon*delta)*phi
          - (1/kon)*2*kp*(koff + kon*(LT + TT)
             - sqrt(-4*kon**2*LT*TT + (koff + kon*(LT + TT))**2))
             *(-1 + alpha)*alpha**4
             *((Gamma_poff + Gamma_pon + YT*delta)*lmbd + (Gamma_yoff + Gamma_yon)*mu)*phi
          + 4*((Gamma_poff + Gamma_pon)*(Gamma_yoff + Gamma_yon) + YT*Gamma_yon*delta)*phi**2
        )

# 3) Simplificar con Sympy
expr_simpl = simplify(expr)
print("\nExpresión simplificada:\n", expr_simpl)

# 4) Sustitución adicional: g1 = (Gamma_yoff + Gamma_yon)
GammaY = symbols("GammaY", positive=True)  # Definimos g1 simbólico
expr_sub = expr_simpl.subs(Gamma_yoff + Gamma_yon, GammaY)
GammaP = symbols("GammaP", positive=True)  # Definimos g2 simbólico
expr_sub1 = expr_sub.subs(Gamma_poff + Gamma_pon, GammaP)

expr_sub_simpl = simplify(expr_sub)
print("\nExpresión tras sustituir g1 = (Gamma_yoff + Gamma_yon) y volver a simplificar:\n", expr_sub_simpl)

# 3) Derivar con respecto a LT
dexpr_dLT = diff(expr_simpl, LT)

# 4) Resolver dexpr_dLT == 0 con respecto a LT
soluciones_LT = solve(sympy.Eq(dexpr_dLT, 0), LT)

#expr_factor = factor(expr_simpl)
#expr_cancel = cancel(expr_simpl)
#expr_rad = radsimp(expr_simpl)

# 4) Convertir a LaTeX
latex_expr = latex(soluciones_LT)

'''
print("Expresión simplificada en LaTeX:")
print(latex_expr)

# 4) Convertir a LaTeX
latex_expr = latex(expr_sub1)

print("Expresión simplificada 1 en LaTeX:")
print(latex_expr)
'''