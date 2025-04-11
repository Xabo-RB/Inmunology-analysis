using Symbolics

@variables X C a1 a3 a4 a6 a7 a8

# Definir la función Tp
Tp = (a1*X - a3*X*C + a4*C^2 - a8*C) / (a6*X - a7*C)

# Derivada parcial con respecto a X
dTp_dX = Symbolics.derivative(Tp, X)

# Derivada parcial con respecto a C
dTp_dC = Symbolics.derivative(Tp, C)

# Mostrar resultados
println("∂Tp/∂X = ", dTp_dX)
println("∂Tp/∂C = ", dTp_dC)

using Symbolics

@variables X C a1 a3 a4 a6 a7 a8

# Definir la función Tp
Tp = (a1*X - a3*X*C + a4*C^2 - a8*C) / (a6*X - a7*C)

# Calcular la Hessiana (segunda derivada parcial)
H = Symbolics.hessian(Tp, [X, C])

# Mostrar la matriz Hessiana
println("Hessiana de Tp:")
display(H)

