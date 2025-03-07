%% KPC
% initial values
x0 = [100, 2, 0, 0, 0, 0, 0]; 
% step size and time interval in days
tspan = 0.0:0.05:600;
% k1 = p[1] = kon,  k3 = p[2], kmenos1 = p[3], w = p[4], k2 = p[5], kmenos2 = p[6]
p = [10, 1, 0.1, 1, 1, 10];

KPC = @(t,y)ODEKPC(t, y, p);
options = odeset('RelTol',1e-5,'AbsTol',1e-6, 'Refine', 1);
[t,x] = ode45(KPC, tspan, x0, options);

plot(tspan,x(:,7));


function dx = ODEKPC(t, x, p)
    dx = zeros(7,1);
    dx(1) = -p(1) * x(1) * x(2) + p(2) * x(4) + p(3) * (x(5) + x(6) + x(7));
    dx(2) = -p(1) * x(1) * x(2) + p(4) * x(7) + p(3) * (x(5) + x(6) + x(7)) ...
            - p(5) * x(3) * x(2) + p(6) * x(4) + p(2) * x(4);
    dx(3) = p(4) * x(7) - p(5) * x(3) * x(2) + p(6) * x(4);
    dx(4) = p(5) * x(3) * x(2) - (p(6) + p(2)) * x(4);
    dx(5) = p(1) * x(1) * x(2) - (p(3) + p(4)) * x(5);
    dx(6) = -p(3) * x(6) - p(4) * x(6) + p(4) * x(5);
    dx(7) = -p(3) * x(7) - p(4) * x(7) + p(4) * x(6);
end