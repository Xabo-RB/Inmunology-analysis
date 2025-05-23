clear
clc

%% Valores parámetros
k1 = 1e-5;
kmenos1 = 5e-2;
w = 9e-2;
k3 = 4e-2;
k2 = 1e-1;
kmenos2 = 5e-2;

%k1 = 10;k3 = 1; kmenos1 = 0.1; w = 1; k2 = 1; kmenos2 = 1;

TT = 3e4;
N = 2;

% Cálculo de α
alpha = w / (kmenos1 + w);
% Cálculo de θ
theta = w * alpha^N * ((1 - alpha) / (1 - alpha^(N + 1)));

% Resultado Emax
EmaxFaro = (alpha + sqrt(alpha^2 + 4 * theta^2 * k3^2 * k2 * (kmenos2 + k3) * TT)) / (2 * theta * k2 * k3);

%% Mis resultados:
EmaxXabo = ((kmenos2+k3)/k2)*TT;
PT = 1000;
EC50Xabo = (2 * k2 * k3 * PT * ((kmenos1 + w)^2 - w * (kmenos1 + 2*w) * (w / (kmenos1 + w))^N)) / ...
       (k1 * w * (w / (kmenos1 + w))^N * (-2 * k2 * k3 * PT + 2 * k2 * kmenos1 * (TT + 2) + kmenos1 * (TT + 2) * (k3 + kmenos2)) + ...
        2 * k2 * k3 * k1 * PT * (kmenos1 + w));

%% Simulación

% initial values
x0_original = [3e4, 5e4, 0, PT, 0, 0, 0, 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.05:600;

p = [1e-5, 4e-2, 5e-2, 9e-2,  1e-1, 5e-2];

% Vector logarítmico de valores de x0(2)
NN = 50; 
x0_values = logspace(0, 7, NN);

max_Tp_values = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPZU(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_Tp_values(i) = max(x(:,3));
    %max_x7_values(i) = x(end,7);
end

% Graficar los resultados
figure;
semilogx(x0_values, max_Tp_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');

Emax = max_Tp_values(end);
Emaxmitad =  Emax/2;
EC50Real = interp1(max_Tp_values, x0_values, Emaxmitad, 'spline');






function dx = ODEKPZU(t, x, p)
    dx = zeros(8,1);
    % k1 = p(1); k3 = p(2); kmenos1 = p(3); w = p(4); k2 = p(5); kmenos2 =
    % p(6)
    dx(1) = -p(1) * x(1) * x(2) + p(2) * x(5) + p(3) * (x(6) + x(7) + x(8));     % T'(t)
    dx(2) = -p(1) * x(1) * x(2) + p(4) * x(8) + p(3) * (x(6) + x(7) + x(8));      % P'(t)
    dx(3) = p(4) * x(8) - p(5) * x(3) * x(4) + p(6) * x(5);                       % Tp'(t)
    dx(4) = -p(5) * x(3) * x(4) + p(6) * x(5) + p(2) * x(5);                     % Q'(t)
    dx(5) = p(5) * x(3) * x(4) - (p(6) + p(2)) * x(5);                           % D'(t)
    dx(6) = p(1) * x(1) * x(2) - (p(3) + p(4)) * x(6);                            % C0'(t)
    dx(7) = -p(3) * x(7) - p(4) * x(7) + p(4) * x(6);                              % C1'(t)
    dx(8) = -p(3) * x(8) - p(4) * x(8) + p(4) * x(7);                              % C2'(t)
end