clear
clc

%% Simulación

TT = 3e4;
LT = 5e4;
N = 5;
% initial values
x0_original = [TT, 5e4, 0, 0, 0, 0, 0, 0, 0, 0, 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.01:80;

%kon = p[1], koff = p[2], kp = p[3],  phi = p[4],   gammaPos = p[5],
%lambda = p[6],  delta = p[7],   YT = p[8],  PT = p[9],  mu = p[10]; gammaNeg = p[11]
p = [1e-5, 0.05, 0.04, 0.04, 1, 100, 100, 100, 2.5, 500, 500];
kon = p(1); koff = p(2); kp = p(3);  phi = p(4);   gammaPos = p(5); lambda = p(6);  delta = p(7);   YT = p(8);  PT = p(9);  mu = p(10); gammaNeg = p(11);

KPC = @(t,y)ODELimIFF1(t, y, p);
[t,x] = ode45(KPC, tspan, x0_original, options);


figure
plot(t, x(:,11), 'LineWidth', 2, 'Color', [0 0 0.6])
hold on
[ymax, idx_max] = max(x(:,11));
t_max = t(idx_max);
h_max = plot(t_max, ymax, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on
box on
xlabel('Time (s)', 'FontSize', 12)
ylabel('P(t)', 'FontSize', 12)
set(gca, 'FontSize', 11)
set(gcf, 'Color', 'w')
legend(h_max, 'Max. response', 'Location', 'best')
% plot(t,x(:,6))

% Vector logarítmico de valores de x0(2)
NN = 50; 
XT_values = logspace(0, 7, NN);

max_CN_values = zeros(size(XT_values));
CN_SS = zeros(size(XT_values));
resultadoCNteorico = zeros(size(XT_values));

psi = kp/(koff/kp);
%alpha = psi;
alpha = (koff / (koff + phi)) * psi^N;

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = XT_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODELimIFF1(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_CN_values(i) = max(x(:,11));
    CN_SS(i) = x(end,11);
    %max_x7_values(i) = x(end,7);

    CT_hat = x(end,3) + x(end,4) + x(end,5) + x(end,6) + x(end,7) + x(end,8) + x(end,9);

    numerador = (gammaPos * (gammaNeg + gammaPos + lambda * alpha * CT_hat) + ...
             delta * (gammaPos + lambda * alpha * CT_hat) * YT) * PT;
    denominador = (gammaNeg + gammaPos + mu * alpha * CT_hat) * ...
              (gammaNeg + gammaPos + lambda * alpha * CT_hat) + ...
              delta * (gammaPos + lambda * alpha * CT_hat) * YT;
    resultadoCNteorico(i) = numerador / denominador;
end

figure
semilogx(XT_values, max_CN_values, '-o')
xlabel('Total ligands')
ylabel('Maximal response')
hold on

figure
semilogx(XT_values, resultadoCNteorico, '-o')
xlabel('Total ligands')
ylabel('Maximal response')
hold on

figure
semilogx(XT_values, CN_SS, '-o')
xlabel('Total ligands')
ylabel('Maximal response')
hold on

figure;
% Primera curva (at Steady-state)
h1 = semilogx(XT_values, CN_SS, '-.', ...
    'LineWidth', 1.2, ...
    'Color', [0.4660 0.6740 0.1880], ...
    'MarkerSize', 1, ...
    'DisplayName', 'at Steady-state');
hold on;
% Segunda curva (Maximum)
h2 = semilogx(XT_values, max_CN_values, '-', ...
    'LineWidth', 1.2, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', 'Maximum');
xlabel('$L_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('$\widehat R\ \mathrm{(response)}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
box off;
legend([h1, h2], 'Interpreter', 'none', 'Location', 'best');


save('resultadosCN_IFF.mat', 'XT_values', 'CN_SS');
save('resultadosCN_IFF1.mat', 'XT_values', 'max_CN_values');
save('resultadosCN_IFF2.mat', 'XT_values', 'resultadoCNteorico');

%%
eMAX_TEORICO = zeros(size(XT_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    LT = XT_values(i);
    
    g1 = (koff + kon * (LT + TT)) ...
         - sqrt( -4 * LT * TT * kon^2 + (koff + kon * (LT + TT))^2 );
    
    Num = 4 * PT * kon^2 * (koff + phi) * ...
        (gammaPos * koff * (gammaNeg + gammaPos) + ...
         gammaPos * phi * (gammaNeg + gammaPos) * gammaPos * YT * delta * koff + ...
         gammaPos * YT * delta * phi) ...
        - 2 * PT * kon * (koff + phi) * psi^4 * kp * lambda * ...
        (gammaPos + YT * delta) * (psi - 1) * g1;
    
    Den = psi^8 * kp^2 * lambda * mu * (psi - 1)^2 * g1^2 + ...
        4 * kon^2 * (gammaPos * YT * delta + (gammaNeg + gammaPos) * (gammaNeg + gammaPos)) * ...
        (koff^2 + 2 * koff * phi + phi^2) - ...
        2 * psi^4 * kon * kp * (psi - 1) * (koff + phi) * ...
        (lambda * (gammaNeg + gammaPos + YT * delta) + mu * (gammaNeg + gammaPos)) * g1;

    eMAX_TEORICO(i) = Num/Den;
end


figure;
semilogx(XT_values, CN_SS, 'o', ...
    'MarkerSize', 4, ...                % Tamaño de marcadores un poco mayor
    'MarkerEdgeColor', 'k');
hold on;
semilogx(XT_values, CN_SS, '-', ...
    'LineWidth', 1, ...
    'Color', [0 0.4470 0.7410]);        % Mismo color azul para la línea
xlabel('$L_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('$\widehat R$ (response at steady-state)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
box off;
% --- Líneas horizontales E_max y E_max/2 ---
[~, idx_max] = max(CN_SS);
maxVal = CN_SS(idx_max);
half_val = maxVal/2;
hLineH = yline(maxVal, ...
    'Color', [0.8500 0.3250 0.0980], ... % naranja rojizo
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$');
% hLineH1 = yline(half_val, ...
%     'Color', [0.4660 0.6740 0.1880], ... % verde
%     'LineStyle', '-', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', '$E_{max}$/2');
x_end = 1e6;
x_start = min(XT_values);  % o lo que tú quieras
hLineH1 = line([x_start x_end], [half_val half_val], ...
    'Color', [0.4660 0.6740 0.1880], ...
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$E_{max}$/2');
%x_half = interp1(CN_SS, XT_values, half_val, 'spline');
mitadInf = CN_SS(1:idx_max);
XT_inf   = XT_values(1:idx_max);
mitadSup = CN_SS(idx_max:end);
XT_sup   = XT_values(idx_max:end);
x_half = interp1(mitadInf, XT_inf, half_val, 'spline');
x_half2 = interp1(mitadSup, XT_sup, half_val, 'spline');
hLineV = xline(x_half, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
hLineV1 = xline(x_half2, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
ylim([0 1.05 * maxVal]); 
legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
hold off;

%% Mi resultado

kon = 5e-5;
koff = 0.01;
kp = 1;
phi = 0.09;
gama = 500;
gammaNeg = 500; gammaPos = 1;
lambda = 0.5;
delta = 50;
YT = 100;
PT = 100;
mu = 2.5;

psi = kp/(kp+koff);
nu = (1 - psi^2)*koff/(koff+phi);

EmaxXabo = PT * (2 + gama + ((lambda*TT/gammaPos) + (delta*YT/gammaPos) * (lambda*TT/gammaPos)) * nu) / ...
         ((lambda*TT/gammaPos) * (mu*TT/gammaPos) * nu^2 + ((lambda*TT/gammaPos) + (lambda*TT/gammaPos) * gama + (mu*TT/gammaPos) + ...
         gama * (mu*TT/gammaPos) + (delta*YT/gammaPos) * (lambda*TT/gammaPos)) * nu + gama^2 + ...
         2 * gama + (delta*YT/gammaPos) + 1);




function dx = ODELimIFF(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(7, 1);

    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5));  % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5));  % T
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);  % C0
    dx(4) = p(3) * x(3) - (p(2) + p(4)) * x(4);  % C1
    dx(5) = p(4) * x(4) - p(2) * x(5);   % C2
    dx(6) = p(5) * (p(8) - x(6)) - p(11) * x(6) + p(6) * x(4) * (p(8) - x(6));  % Y
    dx(7) = p(5) * (p(9) - x(7)) - p(11) * x(7) + p(7) * x(6) * (p(9) - x(7)) - p(10) * x(4) * x(7);    % P
end

function dx = ODELimIFF1(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(7, 1);

    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9));  % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9));  % T
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);  % C0
    dx(4) = p(3) * x(3) - (p(2) + p(3)) * x(4);  % C1
    dx(5) = p(3) * x(4) - (p(2) + p(3)) * x(5);  % C2
    dx(6) = p(3) * x(5) - (p(2) + p(3)) * x(6);  % C3
    dx(7) = p(3) * x(6) - (p(2) + p(3)) * x(7);  % C4
    dx(8) = p(3) * x(7) - (p(2) + p(4)) * x(8);  % C5
    dx(9) = p(4) * x(8) - p(2) * x(9);   % C6
    dx(10) = p(5) * (p(8) - x(10)) - p(11) * x(10) + p(6) * x(8) * (p(8) - x(10));  % Y
    dx(11) = p(5) * (p(9) - x(11)) - p(11) * x(11) + p(7) * x(10) * (p(9) - x(11)) - p(10) * x(8) * x(11);    % P
end