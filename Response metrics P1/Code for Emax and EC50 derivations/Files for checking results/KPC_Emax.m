%% KPC
clear
clc
tic

N = 50; 
% initial values
TT = 3e4; XT = 5e4;

x0_original = [TT, XT, 0, 0, 0, 0, 0];

tolerancia = 1e-4;

% step size and time interval in days
tspan = 0.0:0.5:2000;

kon = 0.00001;  koff = 0.05;    w = 0.09;   k3 = 0.04;  k2 = 0.1;   kmenos2 = 0.05; 

% k1 = p[1] = kon,  k3 = p[2], kmenos1 = p[3], w = p[4], k2 = p[5], kmenos2 = p[6]
%p = [10, 1, 0.1, 1, 1, 10];
%p = [1e-5, 4e-2, 5e-2, 9e-2, 1e-1, 5e-2];
p = [kon, k3, koff, w, k2, kmenos2];

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

figure
plot(tspan,x(:,2));
title('x')
hold on

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
CTss = CT(end);

CT = x(:,5)+x(:,6)+x(:,7);
figure
plot(tspan, CT);
title('CT')
hold on

%Tp = (koff*CTss - kon*CTss^2 + kon*TT*CTss + kon*CTss*XT - kon*TT*XT) / ((k_2*k_3/(k_menos2 + k_3) + kon)*(CTss - XT));



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

% EQN 79 de Tp del overleaf

Num = kon*XT*TT - kon*TT*(1+f1).*CT - (1+f1)*kon*XT.*CT + kon*((1+f1)^2).*CT.^2 - (k3*f1 + koff).*CT;
Den = kon.*(XT - (1+f1).*CT);
TpTeorico = Num./Den;

if (TpTeorico(end) - x(end,3)) < tolerancia
    disp('Tp teórico eqn 3 coincide')
end


% numerador = f1*(k_3 + k_menos2)*(f1*(k_3 + k_menos2) + k_2*(f1*k_3 + k_menos1) - ...
%     (1 + f1) * k_2 * TT + sqrt( (k_2 * (k_menos1 + TT) + f1 * (k_3 + k_2 * k_3 + k_menos2 + k_2 * TT))^2 ));
% 
% denominador = 2 * (1 + f1) * k_2 * (f1 * (k_3 + k_2 * k_3 + k_menos2) + k_2 * k_menos1);
% 
% TpTeoricoMathematica1 = - (f1 * (k_3 + k_menos2) * ...
%        (f1 * (k_3 + k_menos2) + k_2 * (f1 * k_3 + k_menos1) - (1 + f1) * k_2 * TT + ...
%         sqrt((k_2 * (k_menos1 + TT) + f1 * (k_3 + k_2 * k_3 + k_menos2 + k_2 * TT))^2))) ...
%      / (2 * (1 + f1) * k_2 * (f1 * (k_3 + k_2 * k_3 + k_menos2) + k_2 * k_menos1));
% 
% TpTeoricoMathematica2 = (f1 * (k_3 + k_menos2) * ...
%       (-f1 * (k_3 + k_menos2) - k_2 * (f1 * k_3 + k_menos1) + ...
%        (1 + f1) * k_2 * TT + ...
%        sqrt((k_2 * (k_menos1 + TT) + f1 * (k_3 + k_2 * k_3 + k_menos2 + k_2 * TT))^2))) ...
%      / (2 * (1 + f1) * k_2 * (f1 * (k_3 + k_2 * k_3 + k_menos2) + k_2 * k_menos1));
% 
% 
% TpTeoricoMathematica = numerador / denominador;



%% Ecuación haciendo el límite cuando XT tiende a infinito sobre Tp

if x0_original(1) == XT
    TTcopia = TT; 
    TT = XT;
end

% NO FUNCIONA PORQUE ES EL CTss que corresponde al XTi, hay que hallar el
% valor de CT que da el máximo y que dependa de los parámetros. 

CTSS = CT(end);
Num1 = kon * TT * (1 + f1) * CTSS - kon * ((1 + f1)^2) * (CTSS^2) - (k3 * f1 + koff) * CTSS;
Den1 = kon * ((1 + f1) * CTSS);
TpTeoricoXT0 = Num1/Den1;
maxTP = max(x(:,3));

tolerancia1 = 0.5;
if abs((TpTeoricoXT0 - maxTP)) < tolerancia1
    disp('Tp teórico cuando XT=0 coincide con el máximo')
end

Num2 = kon*TT - kon*TT*(1+f1)*CTSS - (1+f1)*kon*CTSS + kon*((1+f1)^2)*CTSS^2 - (k3*f1 + koff)*CTSS;
Den2 = kon*(1 - (1+f1)*CTSS);
TpTeoricoXT1 = Num2/Den2;

if abs((TpTeoricoXT1 - maxTP)) < tolerancia1
    disp('Tp teórico cuando XT=1 coincide con el máximo')
end

figure
plot(tspan, x(:,3), 'DisplayName', 'Tp');
hold on
plot(tspan, TpTeorico, 'DisplayName', 'Tp_teorico');
hold off
legend show
grid on

if x0_original(1) == XT
    TT = TTcopia;
end

%% Usando la ecuación de Faro para Tp

% TpFaro = ((k_menos2+k_3)*w*psi^NN)/(-k_2*( w*psi^NN + ((1-psi^(NN+1))/(1-psi)) ));
% 
% TpFaro1 = -((psi-1) * psi * w * (k_3 + k_menos2)) / (k_2 * (k_3 * (psi^(NN+1) - 1) + (psi-1) * psi * w));
% 
% TpFaro2 = -(sqrt((k_2 * TT * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w) + ...
%     (psi - 1) * psi * w * (k_3 + k_menos2))^2) - ...
%     k_2 * TT * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w) + ...
%     (psi - 1) * psi * w * (k_3 + k_menos2)) / ...
%     (2 * k_2 * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w)); %Primera solución (da TT)
% 
% TpFaro2 = (sqrt((k_2 * TT * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w) + ...
%     (psi - 1) * psi * w * (k_3 + k_menos2))^2) + ...
%     k_2 * TT * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w) - ...
%     (psi - 1) * psi * w * (k_3 + k_menos2)) / ...
%     (2 * k_2 * (k_3 * (psi^(NN+1) - 1) + (psi - 1) * psi * w)); %Segunda solución
% 
% disp(['El Tp de Faro: ', num2str(TpFaro)]);
% disp(['El Tp de Faro: ', num2str(TpFaro1)]);
% disp(['El Tp de Faro: ', num2str(TpFaro2)]);
% disp(['El Tp es: ', num2str(maxTP)]);


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
TpTeoricoXT = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:N
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0

    % Resolver la ODE
    KPC = @(t,y)ODEKPC(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    [max_Tp_values(i), indMax] = max(x(:,3));

    CT = x(:,5)+x(:,6)+x(:,7);
    CTss = CT(indMax);
    CTss = CT(end);

    Num = kon*x0_values(i)*TT - kon*TT*(1+f1)*CTss - (1+f1)*kon*x0_values(i)*CTss + kon*((1+f1)^2)*CTss^2 - (k3*f1 + koff)*CTss;
    Den = kon.*(x0_values(i) - (1+f1)*CTss);
    TpTeoricoXT(i) = Num/Den;
    
    max_Tp_valuesSteadyState(i) = x(end,3);

    max_D_values(i) = x(end,4);

end

%% Tp con punto rojo en el máximo

% Graficar los resultados
% figure;
% semilogx(x0_values, max_Tp_values, '-o');
% xlabel('Total ligands');
% ylabel('Maximal response');
% %title(['Emax para k2 = ', num2str(p(5))]);
% hold on

% figure;
% semilogx(x0_values, max_Tp_values, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0 0 0.6]);
% xlabel('Total ligands','FontSize',12,'FontName','Helvetica');
% ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
% grid on;
% set(gca, 'FontSize',12, 'FontName','Helvetica');
% box off; % estética más limpia sin marco

figure;
semilogx(x0_values, max_Tp_values, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0 0 0.6]);
xlabel('$X_{T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
grid on;
set(gca, 'FontSize',12, 'FontName','Helvetica');
box off;
hold on;
[~, idx_max] = max(max_Tp_values);
semilogx(x0_values(idx_max), max_Tp_values(idx_max), 'ro', 'MarkerSize',9, 'LineWidth',1.5);
hold off;


figure;
semilogx(x0_values, max_Tp_values, 'o', 'MarkerSize', 4, 'Color',[0 0 0.6]);
%semilogx(x0_values, max_Tp_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%semilogx(x0_values, max_Tp_values, 'LineWidth', 1.5, 'Color',[0 0 0.6]);
xlabel('$X_T$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
ylabel('Maximal response','FontSize',12,'FontName','Helvetica');
%grid on;
set(gca, 'FontSize',12, 'FontName','Helvetica');
box off;
hold on;
%TpTeoricoXT
semilogx(x0_values, max_Tp_values, '-', 'LineWidth', 1);
[~, idx_max] = max(max_Tp_values);
half_val = max_Tp_values(idx_max)/2;
hLineH = yline(max_Tp_values(idx_max), ...
    'Color', 'r', ...
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$');
hLineH1 = yline(max_Tp_values(idx_max)/2, ...
    'Color', 'g', ...
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$/2');

x_half = interp1(max_Tp_values, x0_values, half_val, 'spline');
hLineV = xline(x_half, ...
    'Color', [1 0.65 0], ...  % naranja
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
maxVal = max_Tp_values(idx_max);
ylim([0 1.05 * maxVal]);  % Por ejemplo, un 10% extra por encima
legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
hold off;

figure;
semilogx(x0_values, max_Tp_values, 'o', ...
    'MarkerSize', 4, ...                % Tamaño de marcadores un poco mayor
    'MarkerEdgeColor', 'k');
hold on;
semilogx(x0_values, max_Tp_values, '-', ...
    'LineWidth', 1, ...
    'Color', [0 0.4470 0.7410]);        % Mismo color azul para la línea
xlabel('$X_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('$\widehat R$ (response at steady-state)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
box off;
% --- Líneas horizontales E_max y E_max/2 ---
[~, idx_max] = max(max_Tp_values);
maxVal = max_Tp_values(idx_max);
half_val = maxVal/2;
hLineH = yline(maxVal, ...
    'Color', [0.8500 0.3250 0.0980], ... % naranja rojizo
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$');
hLineH1 = yline(half_val, ...
    'Color', [0.4660 0.6740 0.1880], ... % verde
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$/2');
x_half = interp1(max_Tp_values, x0_values, half_val, 'spline');
hLineV = xline(x_half, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');

ylim([0 1.05 * maxVal]); 
legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
hold off;

% figure;
% semilogx(x0_values, TpTeoricoXT, '-', 'LineWidth', 1.5);
% set(gca, 'FontSize',12, 'FontName','Helvetica');
% hold on;
% 
% figure;
% semilogx(x0_values, max_Tp_valuesSteadyState, '-', 'LineWidth', 1.5);
% set(gca, 'FontSize',12, 'FontName','Helvetica');
% hold on;
% semilogx(x0_values, max_Tp_values, '-', 'LineWidth', 1.5);

CTparaTpmax = max_CT_values(idx_max);
disp(['El valor de CT cuando Tp es max:', CTparaTpmax])

%% Doble gráfica

figure;
% Primer eje (izquierda) - primera variable
yyaxis left
semilogx(x0_values, max_Tp_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 0.6]);
xlabel('$X_{T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
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
ylabel('$\widehat{C_T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
%title('$\widehat{C_T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
hold on

figure;
semilogx(x0_values, max_CT_values, '-o', 'LineWidth',1.5,'MarkerSize',6,'Color',[0 0 0.6]);
xlabel('$X_{T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
ylabel('$\widehat{C_T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
set(gca,'FontSize',12,'FontName','Helvetica');
grid on;
box off;
% Destacar máximo claramente:
hold on;
[~, idx_max] = max(max_CT_values);
semilogx(x0_values(idx_max), max_CT_values(idx_max), 'ro','MarkerSize',9,'LineWidth',1.5);
yline(3e4, '--', '$T_T$', 'Interpreter','latex', 'FontSize',13, ...
    'LineWidth', 1.5, 'LabelVerticalAlignment','bottom','Color',[0 0.5 0]);
ylim([min(max_CT_values)*0.9, 3.3e4]);
hold off;


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

%% Tiempo y Xt

% step size and time interval in days
tfinal = 50:50:3000;

% Vector logarítmico de valores de x0(2)
x0_values = logspace(0, 7, N);
XtMax = zeros(size(tfinal));
max_Tp_values = zeros(size(x0_values));

tspan = 0.0:0.05:300;

for j = 1:length(tfinal)
    
    tspan = 0.0:0.05:tfinal(j);
    % Bucle sobre cada valor de x0(2)
    for i = 1:N
        x0 = x0_original;
        x0(2) = x0_values(i); % Modificamos el segundo valor de x0
        
        % Resolver la ODE
        KPC = @(t,y)ODEKPC(t, y, p);
        [t, x] = ode23s(KPC, tspan, x0, options);
        
        max_Tp_values(i) = max(x(:,3));
        %XtMax(i) = x(indiceMax,2);
    
    end
    
    [a , indiceMax] = max(max_Tp_values);
    XtMax(j) = x0_values(indiceMax);
end

figure
plot(tfinal, XtMax,'-o');
grid on
hold on

figure;
plot(tfinal, XtMax, '-o', 'LineWidth',1.5,'MarkerSize',5,'Color',[0 0 0.6]);
xlabel('Time','FontSize',12,'FontName','Helvetica');
ylabel('$X_{T}$','Interpreter','latex','FontSize',14,'FontName','Helvetica');
set(gca,'FontSize',12,'FontName','Helvetica');
grid on;
box off;
hold on;

%yline(0,'--','LineWidth',1.2,'Color',[0.8 0 0]);
% % Ajustar límite inferior claramente por debajo de cero
% ylim([-1, max(XtMax)*1.05]);

% Leyenda ajustada (sin steady-state)
% legend({'Response'},...
%     'FontSize',12,'FontName','Helvetica','Location','northeast','box','off');
% hold off;


% Guardar los resultados
save('resultadosXt.mat', 'XtMax', 'tfinal');
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



