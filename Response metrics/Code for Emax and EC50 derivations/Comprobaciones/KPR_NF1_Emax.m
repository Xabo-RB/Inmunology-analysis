clear
clc

TT = 3e4;
ST = 6e5;
x0 = [TT, 5e4, 0, 0, 0, 0, 0, 0, ST]; 
%x0 = [TT, 5e4, 0, 0, 0, ST]; 

% step size and time interval in days
d = 1.0e-16; 
tspan = 0.0:0.05:200;

kon = 1 / 1e5;
koff = 0.03;
kp = 0.13;
b = 0.04;
gama = 1 / 1e6;
alpha = 1 / (5 * 1e2);
beta = 1;


% kon koff kp gama b beta alpha ST
p = [kon koff kp gama b beta alpha ST];

%p = [5e-5, 0.01, 1, 4.4e-4, 0.04, 1, 2e-4, 6e5];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);
Neg1 = @(t,y)ODEKPRNegFeed(t, y, p);
[t,x] = ode45(Neg1, tspan, x0, options);

% figure
% plot(tspan,x(:,5))
% hold on
% 
% figure
% plot(tspan,x(:,6))
% hold on

figure
plot(tspan,x(:,5))
hold on


%P_hat = x(end,6);
P = x(:,6); 
C2 = x(:,5);
k1 = koff + b + gama.*P;
psi = kp ./ (kp + k1);
CT = ( (1./psi.^2)- (b + gama.*P)./kp).*C2;

P = x(end,6); 
C2 = x(end,5);
k1 = koff + b + gama*P;
psi = kp / (kp + k1);
CThat = ( (1./psi^2)- (b + gama*P)./kp)*C2;
CThat_real=x(end,4)+x(end,3)+x(end,5);

%% Simulación
NN = 50; 
x0_values = logspace(0, 7, NN);
x0_original = x0;

max_CN_values = zeros(size(x0_values));
max_CN_SS = zeros(size(x0_values));
P_SS = zeros(size(x0_values));
respuestaNF2 = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i);
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPRNegFeed(t, y, p);
    [t, x] = ode45(KPC, tspan, x0, options);
    
    %   N = 5
    max_CN_values(i) = max(x(:,8));
    max_CN_SS(i) = x(end,8);
    P_SS(i) = x(end,9);

    respuestaNF2(i) = max(x(:,8)) + max(x(:,7))*exp(-2) + max(x(:,6))*exp(-4);

    %   N = 2
    % max_CN_values(i) = max(x(:,5));
    % max_CN_SS(i) = x(end,5);
    % P_SS(i) = x(end,6);
    % 
    % respuestaNF2(i) = max(x(:,5)) + max(x(:,4))*exp(-2) + max(x(:,3))*exp(-4);
end

figure;
semilogx(x0_values, max_CN_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');
hold on
semilogx(x0_values, max_CN_SS, '-o');
% 
figure;
semilogx(x0_values, respuestaNF2, '-o');
xlabel('Total ligands');
ylabel('Maximal response');
% 
EmaxReal = max_CN_SS(end);

%%      N = 5
[~, idx_max] = max(respuestaNF2);
maxVal = respuestaNF2(idx_max);
half_val = maxVal/2;

mitadInf = respuestaNF2(1:idx_max);
mitadInfLT = x0_values(1:idx_max);

mitadSup = respuestaNF2(21:25);
mitadSupLT = x0_values(21:25);

mitadSup1 = respuestaNF2(27:32);
mitadSupLT1 = x0_values(27:32);

x_half1 = interp1(mitadInf, mitadInfLT, half_val, 'spline');
x_half2 = interp1(mitadSup, mitadSupLT, half_val, 'spline');
x_half3 = interp1(mitadSup1, mitadSupLT1, half_val, 'spline');
%x_half2 = x0_values(27);

figure;
semilogx(x0_values, respuestaNF2, 'o', ...
    'MarkerSize', 4, ...                % Tamaño de marcadores un poco mayor
    'MarkerEdgeColor', 'k');
hold on;
semilogx(x0_values, respuestaNF2, '-', ...
    'LineWidth', 1, ...
    'Color', [0 0.4470 0.7410]);        % Mismo color azul para la línea
xlabel('$L_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('$\widehat R$ (response at steady-state)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
box off;
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
hLineV = xline(x_half1, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
hLineV1 = xline(x_half2, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
hLineV1 = xline(x_half3, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
ylim([0 1.05 * maxVal]); 
legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
hold off;

%% N = 2
% figure;
% semilogx(x0_values, respuestaNF2, 'o', ...
%     'MarkerSize', 4, ...                % Tamaño de marcadores un poco mayor
%     'MarkerEdgeColor', 'k');
% hold on;
% semilogx(x0_values, respuestaNF2, '-', ...
%     'LineWidth', 1, ...
%     'Color', [0 0.4470 0.7410]);        % Mismo color azul para la línea
% xlabel('$L_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
% ylabel('$\widehat R$ (response at steady-state)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Helvetica');
% set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
% box off;
% [~, idx_max] = max(respuestaNF2);
% maxVal = respuestaNF2(idx_max);
% half_val = maxVal/2;
% hLineH = yline(maxVal, ...
%     'Color', [0.8500 0.3250 0.0980], ... % naranja rojizo
%     'LineStyle', '-', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', '$E_{max}$');
% hLineH1 = yline(half_val, ...
%     'Color', [0.4660 0.6740 0.1880], ... % verde
%     'LineStyle', '-', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', '$E_{max}$/2');
% x_half = interp1(respuestaNF2, x0_values, half_val, 'spline');
% hLineV = xline(x_half, ...
%     'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
%     'LineStyle', '-', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', '$EC_{50}$');
% ylim([0 1.05 * maxVal]); 
% legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
% hold off;

%% Mis resultados

E_maxXabo = (kp * (koff^2 + koff * kp + kp^2 + b * (koff + kp))) / ...
        (kp * (koff + kp)^2 + b^2 * (koff + 2 * kp) + b * (koff^2 + 4 * koff * kp + 2 * kp^2)) ...
        * TT;

P_hat = 10;
P_hat = P_SS(end);

E_maxXabo1 = (kp * (kp - koff * (psi - 1))) / (kp * (koff + kp) + (b + gama * P_hat) * (koff - kp * (psi - 2))) ...
        * TT;



function dx = ODEKPRNegFeed_N_2(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(6, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5)); % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5)); % T
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(6)) * x(4); % C0
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(6)) * x(4) + (p(5) + p(4) * x(6)) * x(5); % C1
    dx(5) = p(3) * x(4) - (p(2) + p(5) + p(4) * x(6)) * x(5); % C2
    dx(6) = p(7) * x(4) * (p(8) - x(6)) - p(6) * x(6); % P
end


function dx = ODEKPRNegFeed(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(9, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8));
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8));
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(9)) * x(4);
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(4) + (p(5) + p(4) * x(9)) * x(5); % C1
    dx(5) = p(3) * x(4) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(5) + (p(5) + p(4) * x(9)) * x(6); % C2
    dx(6) = p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7); % C3
    dx(7) = p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8); % C4
    dx(8) = p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8); % C5
    dx(9) = p(7) * x(4) * (p(8) - x(9)) - p(6) * x(9);
end