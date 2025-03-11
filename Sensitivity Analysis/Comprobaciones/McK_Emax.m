clear
clc

%% Simulación

TT = 3e4;
% initial values
x0_original = [TT, 5e4 0, 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.05:50;

%kon = p[1], koff = p[2], kp = p[3]
p = [5e-5, 0.01, 1];

% Vector logarítmico de valores de x0(2)
NN = 50; 
x0_values = logspace(0, 7, NN);

max_CN_values = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPRmcK(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_CN_values(i) = max(x(:,4));
    %max_x7_values(i) = x(end,7);
end

figure;
semilogx(x0_values, max_CN_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');

kp = 1; koff = 0.01;
phi = kp/(kp+koff);
Emax = phi^2*TT;
EmaxSim = max_CN_values(end);

function dx = ODEKPRmcK(t, x, p)
    dx = zeros(4,1);
    dx(1) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(2) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);
    dx(4) = p(3) * x(3) - p(2) * x(4);
end