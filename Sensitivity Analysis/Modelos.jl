function ODE(dx, x, p, t) # CARRGO model
    dx[1] = p[4] * x[1] * (1 - x[1] / p[5]) - p[1] * x[1] * x[2]
    dx[2] = p[2]* x[1] * x[2] - p[3] * x[2]
end

function ODEOccupancy(dx, x, p, t) # Occuppancy model
    dx[1] = p[1]*p[2]*p[3] - p[4]*x[1]
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
