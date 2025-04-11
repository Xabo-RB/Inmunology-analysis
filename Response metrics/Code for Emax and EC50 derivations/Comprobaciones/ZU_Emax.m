clear
clc


% k1 = p(1); k3 = p(2); kmenos1 = p(3); w = p(4); k2 = p(5); kmenos2 =
% p(6)

k1 = 1e-5;
kmenos1 = 5e-2;
w = 9e-2;
k3 = 4e-2;
k2 = 1e-1;
kmenos2 = 5e-2;

TT = 3e4;
PT = 10^5;

c = (sqrt(k2^2 * (k3 * PT + w * (PT - TT))^2 - 2 * k2 * w * (k3 + kmenos2) * (k3 * PT + w * (PT + TT)) + w^2 * (k3 + kmenos2)^2) / (2 * k2 * w));
d = - (PT * (k3 + w) / (2 * w)) + (k3 / (2 * k2)) + (kmenos2 / (2 * k2));
EmaxXabo = (TT / 2) + d - c;

alpha = (w * k2 * k3 * TT) - (k2 * k3 * (k3 + w) * PT) - (w * k3 * (kmenos2 + k3));
EmaxFaro = (alpha + sqrt(alpha^2 + 4 * w^2 * k3^2 * k2 * (kmenos2 + k3) * TT)) / (2 * w * k2 * k3);



%% Simulación
p = [1e-5, 4e-2, 5e-2, 9e-2, 1e-1, 5e-2];
LT = 1e4;
x0_original = [3e4, LT, 0, 0, 0, PT];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.05:1000;

% % Vector logarítmico de valores de x0(2)
% NN = 50; 
% x0_values = logspace(0, 7, NN);
% 
% max_Tp_values = zeros(size(x0_values));
% 
% for i = 1:NN
%     x0 = x0_original;
%     x0(2) = x0_values(i); % Modificamos el segundo valor de x0
% 
%     % Resolver la ODE
%     KPC = @(t,y)ODEZU(t, y, p);
%     [t, x] = ode23s(KPC, tspan, x0, options);
% 
%     max_Tp_values(i) = max(x(:,5));
%     %max_x7_values(i) = x(end,7);
% end
% 
% figure;
% semilogx(x0_values, max_Tp_values, '-o');
% xlabel('Total ligands');
% ylabel('Maximal response');
% hold on

KPC = @(t,y)ODEZU(t, y, p);
[t,x] = ode23s(KPC, tspan, x0_original, options);

figure
plot(tspan,x(:,5));
hold on

figure
plot(tspan,x(:,2));
hold on
plot(tspan,x(:,6));
plot(tspan,LT-x(:,2));
legend('L(t)', 'P(t)', 'LT-L(t)')

function dx = ODEZU(t, x, p)
    dx = zeros(6,1);
    % k1 = p(1); k3 = p(2); kmenos1 = p(3); w = p(4); k2 = p(5); kmenos2 =
    % p(6)
    dx(1) = -p(1) * x(1) * x(2) + p(2) * x(4) + p(3) * x(3);                % T'(t)
    dx(2) = -p(1) * x(1) * x(2) + p(4) * x(3) + p(3) * x(3);                  % L'(t)
    dx(3) = p(1) * x(1) * x(2) - (p(3) + p(4)) * x(3);                        % C0'(t)
    dx(4) = p(5) * x(5) * x(6) - (p(6) + p(2)) * x(4);                       % D'(t)
    dx(5) = -p(5) * x(5) * x(6) + p(6) * x(4) + p(4)* x(3);                  % Tp'(t)
    dx(6) = -p(5) * x(5) * x(6) + p(6) * x(4) + p(2) * x(4);                 % Q'(t)
end
