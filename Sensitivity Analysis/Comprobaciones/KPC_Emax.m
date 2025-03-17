%% KPC
clear
clc
tic

N = 50; 
% initial values
TT = 3e4; XT = 200;
x0_original = [TT, XT, 0, 0, 0, 0, 0];

tolerancia = 1e-4;

% step size and time interval in days
tspan = 0.0:0.05:300;

% k1 = p[1] = kon,  k3 = p[2], kmenos1 = p[3], w = p[4], k2 = p[5], kmenos2 = p[6]
%p = [10, 1, 0.1, 1, 1, 10];
%p = [1e-5, 4e-2, 5e-2, 9e-2, 1e-1, 5e-2];
p = [0.00001, 0.04, 0.05, 0.09, 0.1, 0.05];

kon = p(1);
k_3 = p(2); k3 = p(2);
k_menos1 = p(3);
w = p(4);
k_2 = p(5);
k_menos2 = p(6);
koff = k_menos1;
NN = 2;


% tic
% KPC = @(t,y)ODEKPC(t, y, p);
% options = odeset('RelTol',1e-5,'AbsTol',1e-6, 'Refine', 1);
% [t,x] = ode45(KPC, tspan, x0, options);
% toc

KPC = @(t,y)ODEKPC(t, y, p);
options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);
[t,x] = ode23s(KPC, tspan, x0_original, options);

figure
plot(tspan,x(:,4));
title('D')
hold on

CT = x(:,5)+x(:,6)+x(:,7);
figure
plot(tspan, CT);
title('CT')
hold on
yline(TT, '--r', ['TT = ', num2str(TT)]); % Línea horizontal en TT con etiqueta
CTss = CT(end);

Tp = (koff*CTss - kon*CTss^2 + kon*TT*CTss + kon*CTss*XT - kon*TT*XT) / ((k_2*k_3/(k_menos2 + k_3) + kon)*(CTss - XT));



%% Coincide Tp teórico con el simulado?

TpTeorico1 = TT - x(:,4) - CT - (k3.*x(:,4) + koff.*CT) ./ (kon.*(XT - CT - x(:,4)));

if (TpTeorico1(end) - x(end,3)) < tolerancia
    disp('Tp teórico eqn 1 coincide')
end

psi = w / (w + koff);
f1 = (w/k3) *( (psi^NN - psi^(NN+1)) / (1 - psi^(NN+1)) );

TpTeorico2 = TT - f1.*CT - CT - (k3.*f1.*CT + koff.*CT) ./ (kon.*(XT - CT - f1.*CT));

if (TpTeorico2(end) - x(end,3)) < tolerancia
    disp('Tp teórico eqn 2 coincide')
end


Num = kon*XT*TT - kon*TT*(1+f1).*CT - (1+f1)*kon*XT.*CT + kon*((1+f1)^2).*CT.^2 - (k3*f1 + koff).*CT;
Den = kon.*(XT - (1+f1).*CT);
TpTeorico = Num./Den;

if (TpTeorico(end) - x(end,3)) < tolerancia
    disp('Tp teórico eqn 3 coincide')
end

% figure
% plot(tspan, x(:,3), 'DisplayName', 'Tp');
% hold on
% plot(tspan, TpTeorico, 'DisplayName', 'Tp_teorico');
% hold off
% legend show
% grid on



%% SE CUMPLE QUE D = w/k3 CN y D = ... CT
figure
plot(tspan, x(:,4), 'DisplayName', 'D');
hold on
plot(tspan, x(:,7), 'DisplayName', 'Cn');
plot(tspan, x(:,7)*(w/k3), 'DisplayName', 'DconCN');
plot(tspan, CT*f1, 'DisplayName', 'DconCT');
hold off
legend show
grid on

%%

XTsimulado = x(:,2)+x(:,4)+CT;
TTsimulado = x(:,1) + x(:,3) + x(:,4) + CT;

if all(abs(XTsimulado - XT) < tolerancia) && all(abs(TTsimulado - TT) < tolerancia)
    disp('Correcto XT y TT en todos los elementos')
else
    disp('No es correcto XT y TT en al menos un elemento')
end



%%

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

%%

% Graficar los resultados
figure;
semilogx(x0_values, max_Tp_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');
%title(['Emax para k2 = ', num2str(p(5))]);
hold on

figure;
semilogx(x0_values, max_Tp_values, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0 0 0.6]);
xlabel('Total ligands','FontSize',12,'FontName','Helvetica');
ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
grid on;
set(gca, 'FontSize',12, 'FontName','Helvetica');
box off; % estética más limpia sin marco

figure;
semilogx(x0_values, max_Tp_values, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0 0 0.6]);
xlabel('Total ligands','FontSize',12,'FontName','Helvetica');
ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
grid on;
set(gca, 'FontSize',12, 'FontName','Helvetica');
box off;
hold on;
[~, idx_max] = max(max_Tp_values);
semilogx(x0_values(idx_max), max_Tp_values(idx_max), 'ro', 'MarkerSize',9, 'LineWidth',1.5);
hold off;

%% Doble gráfica
figure;
% Primer eje (izquierda) - primera variable
yyaxis left
semilogx(x0_values, max_Tp_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 0.6]);
xlabel('Total ligands','FontSize',12,'FontName','Helvetica');
ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
set(gca, 'FontSize',12, 'FontName','Helvetica');
box off;
grid on;
hold on;
[~, idx_max] = max(max_Tp_values);
semilogx(x0_values(idx_max), max_Tp_values(idx_max), 'ro', 'MarkerSize',9, 'LineWidth',1.5);

% Segundo eje (derecha) - segunda variable
yyaxis right
semilogx(x0_values, max_D_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0.5 0]);
ylabel('$\widehat{D}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');

% Colorear ambos ejes en negro (sin colores)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

grid on;
box off;

%%

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
%title(['Valor k2 = ', num2str(p(5))]);
hold on

% figure;
% semilogx(x0_values, max_Tp_valuesSteadyState, '-o');
% xlabel('Total ligands');
% ylabel('Maximal response SS');

%%


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



