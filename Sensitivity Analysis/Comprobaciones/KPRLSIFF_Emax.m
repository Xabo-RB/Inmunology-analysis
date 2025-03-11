clear
clc

%% Simulación

TT = 3e4;
% initial values
x0_original = [TT, 5e4, 0, 0, 0, 0, 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.05:100;

%kon = p[1], koff = p[2], kp = p[3],  phi = p[4],   gamma = p[5],   lambda = p[6],  delta = p[7],   YT = p[8],  PT = p[9],  mu = p[10]
p = [5e-5, 0.01, 1, 0.09, 1, 0.5, 50, 100, 100, 2.5, 500];

% Vector logarítmico de valores de x0(2)
NN = 50; 
x0_values = logspace(0, 7, NN);

max_CN_values = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODELimIFF(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_CN_values(i) = max(x(:,7));
    %max_x7_values(i) = x(end,7);
end

figure;
semilogx(x0_values, max_CN_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');

%% Mi resultado

kon = 5e-5;
koff = 0.01;
kp = 1;
phi = 0.09;
gama = 1;
gammaNeg = gama; gammaPos = gama;
lambda = 0.5;
delta = 50;
YT = 100;
PT = 100;
mu = 2.5;

psi = kp/(kp+koff);
nu = (1 - psi^2)*koff/(koff+phi);

EmaxXabo = P_T * (2 + gama + ((lambda*TT/gama) + (delta*YT/gama) * (lambda*TT/gama)) * nu) / ...
         ((lambda*TT/gama) * (mu*TT/gama) * nu^2 + ((lambda*TT/gama) + (lambda*TT/gama) * gama + (mu*TT/gama) + ...
         gama * (mu*TT/gama) + (delta*YT/gama) * (lambda*TT/gama)) * nu + gama^2 + ...
         2 * gama + (delta*YT/gama) + 1);




function dx = ODELimIFF(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(7, 1);

    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5));  % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5));  % T
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);  % C0
    dx(4) = p(3) * x(3) - (p(2) + p(4) * p(3)) * x(4);  % C1
    dx(5) = p(4) * p(3) * x(4) - p(2) * x(5);   % C2
    dx(6) = p(5) * (p(8) - x(6)) - p(11) * x(6) + p(6) * x(4) * (p(8) - x(6));  % Y
    dx(7) = p(5) * (p(9) - x(7)) - p(11) * x(7) + p(7) * x(6) * (p(9) - x(7)) - p(10) * x(4) * x(7);    % P
end