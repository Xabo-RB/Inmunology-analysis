%% KPC
clear
clc
tic

N = 50; 
% initial values
TT = 3e4;
x0_original = [3e4, 5e4, 0, 0, 0, 0, 0];

% step size and time interval in days
tspan = 0.0:0.05:300;

% k1 = p[1] = kon,  k3 = p[2], kmenos1 = p[3], w = p[4], k2 = p[5], kmenos2 = p[6]
%p = [10, 1, 0.1, 1, 1, 10];
%p = [1e-5, 4e-2, 5e-2, 9e-2, 1e-1, 5e-2];
p = [0.00001, 0.04, 0.05, 0.09, 0.1, 0.05];


% tic
% KPC = @(t,y)ODEKPC(t, y, p);
% options = odeset('RelTol',1e-5,'AbsTol',1e-6, 'Refine', 1);
% [t,x] = ode45(KPC, tspan, x0, options);
% toc

KPC = @(t,y)ODEKPC(t, y, p);
options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);
[t,x] = ode23s(KPC, tspan, x0_original, options);

figure
plot(tspan,x(:,3));
title('Tp')
hold on

CT = x(:,5)+x(:,6)+x(:,7);
figure
plot(tspan, CT);
title('CT')
hold on
yline(TT, '--r', ['TT = ', num2str(TT)]); % Línea horizontal en TT con etiqueta

% Vector logarítmico de valores de x0(2)
x0_values = logspace(0, 7, N);
max_Tp_valuesSteadyState = zeros(size(x0_values));
max_Tp_values = zeros(size(x0_values));
max_CT_values = zeros(size(x0_values));
max_D_values = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:N
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPC(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_Tp_values(i) = max(x(:,3));
    max_Tp_valuesSteadyState(i) = x(end,3);
    
    CT = x(:,5)+x(:,6)+x(:,7);
    max_CT_values(i) = max(CT);

    max_D_values(i) = x(end,4);

end

% Graficar los resultados
figure;
semilogx(x0_values, max_Tp_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');
title(['Emax para k2 = ', num2str(p(5))]);
hold on

figure;
semilogx(x0_values, max_CT_values, '-o');
xlabel('Total ligands');
ylabel('CT');
title(['Valor k2 = ', num2str(p(5))]);
hold on

figure;
semilogx(x0_values, max_D_values, '-o');
xlabel('Total ligands');
ylabel('D');
title(['Valor k2 = ', num2str(p(5))]);
hold on

% figure;
% semilogx(x0_values, max_Tp_valuesSteadyState, '-o');
% xlabel('Total ligands');
% ylabel('Maximal response SS');

kon = p(1);
k_3 = p(2);
k_menos1 = p(3);
w = p(4);
k_2 = p(5);
k_menos2 = p(6);
koff = k_menos1;

NN = 2;

T_T = x(1);

psi = w / (w + koff);
resultado = (w * ((psi^NN - psi^(NN+1)) / (1 - psi^(NN+1)))) / (2 * ((k_menos2 * k_2) / (k_menos2 + k_3) - k_2)) + ...
            (T_T / (2 * (1 + (k_2 * k_3) / (kon * (k_menos2 + k_3)))));

disp(['El Emax es: ', num2str(resultado)]);

toc

function dx = ODEKPC(t, x, p)
    dx = zeros(7,1);
    dx(1) = -p(1) * x(1) * x(2) + p(2) * x(4) + p(3) * (x(5) + x(6) + x(7)); % T(t)
    dx(2) = -p(1) * x(1) * x(2) + p(4) * x(7) + p(3) * (x(5) + x(6) + x(7)) ...
            - p(5) * x(3) * x(2) + p(6) * x(4) + p(2) * x(4); % X(t)
    dx(3) = p(4) * x(7) - p(5) * x(3) * x(2) + p(6) * x(4); % Tp
    dx(4) = p(5) * x(3) * x(2) - (p(6) + p(2)) * x(4); % D
    dx(5) = p(1) * x(1) * x(2) - (p(3) + p(4)) * x(5); % C0
    dx(6) = -p(3) * x(6) - p(4) * x(6) + p(4) * x(5); % C1
    dx(7) = -p(3) * x(7) - p(4) * x(7) + p(4) * x(6); % C2
end



