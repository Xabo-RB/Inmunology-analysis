clear
clc

TT = 10^3;

x0 = [TT, 5e4, 0, 0, 0, 0, 0, 0, 0]; 
% step size and time interval in days
d = 1.0e-16; 
tspan = 0.0:0.05:200;

kon = 1 / 1e5;
koff = 0.05;
kp = 0.04;
b = 0.04;
gama = 1 / 1e6;
ST = 6e5;
alpha = 1 / (5 * 1e2);
beta = 1;


% kon koff kp gama b beta alpha ST
p = [kon koff kp gama b beta alpha ST];

%p = [5e-5, 0.01, 1, 4.4e-4, 0.04, 1, 2e-4, 6e5];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);
Neg1 = @(t,y)ODEKPRNegFeed_N_2(t, y, p);
[t,x] = ode45(Neg1, tspan, x0, options);

figure
plot(tspan,x(:,5))
hold on


%% Simulación
NN = 50; 
x0_values = logspace(0, 7, NN);
x0_original = x0;

max_CN_values = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i);
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPRNegFeed_N_2(t, y, p);
    [t, x] = ode45(KPC, tspan, x0, options);
    
    max_CN_values(i) = max(x(:,5));
    %max_x7_values(i) = x(end,7);
end

figure;
semilogx(x0_values, max_CN_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');

EmaxReal = max_CN_values(end);

%% Mis resultados

E_maxXabo = (kp * (koff^2 + koff * kp + kp^2 + b * (koff + kp))) / ...
        (kp * (koff + kp)^2 + b^2 * (koff + 2 * kp) + b * (koff^2 + 4 * koff * kp + 2 * kp^2)) ...
        * TT;

P_hat = 6e5;
k1 = koff + b + gama * P_hat;
psi = kp / (kp + k1);

E_maxXabo1 = (kp * (kp - koff * (psi - 1))) / (kp * (koff + kp) + (b + gama * P_hat) * (koff - kp * (psi - 2))) ...
        * TT;



function dx = ODEKPRNegFeed_N_2(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(9, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4));
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4));
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(6)) * x(4);
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(6)) * x(4) + (p(5) + p(4) * x(6)) * x(5);
    dx(5) = p(3) * x(4) - (p(2) + p(5) + p(4) * x(6)) * x(5);
    dx(6) = p(7) * x(4) * (p(8) - x(6)) - p(6) * x(6);
end


function dx = ODEKPRNegFeed(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(9, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8));
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8));
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(9)) * x(4);
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(4) + (p(5) + p(4) * x(9)) * x(5);
    dx(5) = p(3) * x(4) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(5) + (p(5) + p(4) * x(9)) * x(6);
    dx(6) = p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7);
    dx(7) = p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8);
    dx(8) = p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8);
    dx(9) = p(7) * x(4) * (p(8) - x(9)) - p(6) * x(9);
end