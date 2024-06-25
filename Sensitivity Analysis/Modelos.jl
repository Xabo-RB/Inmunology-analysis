function ODE(dx, x, p, t) # CARRGO model
    dx[1] = p[4] * x[1] * (1 - x[1] / p[5]) - p[1] * x[1] * x[2]
    dx[2] = p[2]* x[1] * x[2] - p[3] * x[2]
end

function ODEOccupancy(dx, x, p, t) # Occuppancy model
    dx[1] = p[1]*x[2]*x[3] - p[2]*x[1]                      #C0
    dx[2] = - p[1] * x[2] * x[3] + p[2]*x[1]                # P
    dx[3] = - p[1] * x[3] * x[3] + p[2]*x[1]                # T

end

function ODEKPRmcK(dx, x, p, t) # McKeithan
    dx[1] = - p[1] * x[1] * x[2] + p[2]*x[3] + p[2]*x[4]
    dx[2] = - p[1] * x[1] * x[2] + p[2]*x[3] + p[2]*x[4]
    dx[3] = p[1] * x[1] * x[2] - (p[2] + p[3])*x[3]
    dx[4] = p[3]*x[3] - p[2]*x[4] #C1 = R

    #   p = complex([5e-5, 0.01, 0.01]);
    #   kon = 5e-5 // koff = 0.01 // kp = 0.01
    #   x0 = complex([100, 2e4, 0, 0]);
    #   P = 100 // T = 20.000 // C0 = 0 // C1 = 0
end

# con esta sí que obtengo el resultado del Paper
function ODEKPRmcK10(dx, x, p, t) # McKeithan
    dx[1] = - p[1] * x[1] * x[2] + p[2]*x[3] + p[2]*x[4] + p[2]*x[5] + p[2]*x[6] + p[2]*x[7] + p[2]*x[8] + p[2]*x[9] + p[2]*x[10] + p[2]*x[11] + p[2]*x[12] + p[2]*x[13] # P
    dx[2] = - p[1] * x[1] * x[2] + p[2]*x[3] + p[2]*x[4] + p[2]*x[5] + p[2]*x[6] + p[2]*x[7] + p[2]*x[8] + p[2]*x[9] + p[2]*x[10] + p[2]*x[11] + p[2]*x[12] + p[2]*x[13] # T
    dx[3] = p[1] * x[1] * x[2] - (p[2] + p[3])*x[3] # C0
    dx[4] = p[3]*x[3] - (p[2] + p[3])*x[4] # C1
    dx[5] = p[3]*x[4] - (p[2] + p[3])*x[5] # C2
    dx[6] = p[3]*x[5] - (p[2] + p[3])*x[6] # C3
    dx[7] = p[3]*x[6] - (p[2] + p[3])*x[7] # C4
    dx[8] = p[3]*x[7] - (p[2] + p[3])*x[8] # C5
    dx[9] = p[3]*x[8] - (p[2] + p[3])*x[9] # C6
    dx[10] = p[3]*x[9] - (p[2] + p[3])*x[10] # C7
    dx[11] = p[3]*x[10] - (p[2] + p[3])*x[11] # C8
    dx[12] = p[3]*x[11] - (p[2] + p[3])*x[12] # C9
    dx[13] = p[3]*x[12] - p[2]*x[13] # C10
    #dx[13] = p[3]*x[12] # C10

end
