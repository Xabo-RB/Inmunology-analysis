clear
clc

%% Simulación

TT = 3e4;
% initial values
x0_original = [TT, 5e4 0, 0, 0, 0, 0, 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Refine', 1);

% step size and time interval in days
tspan = 0.0:0.05:50;

%kon = p[1], koff = p[2], kp = p[3]
kp = 1; koff = 0.01; kon = 5e-5;
p = [5e-5, 0.01, 1];

phi = kp/(kp+koff);
Emax = (phi^5)*TT;

% Vector logarítmico de valores de x0(2)
NN = 50; 
x0_values = logspace(0, 7, NN);

max_CN_values = zeros(size(x0_values));
max_CN_valuesSS = zeros(size(x0_values));
RespuestaSS = zeros(size(x0_values));

% Bucle sobre cada valor de x0(2)
for i = 1:NN
    x0 = x0_original;
    x0(2) = x0_values(i); % Modificamos el segundo valor de x0
    
    % Resolver la ODE
    KPC = @(t,y)ODEKPRmcK(t, y, p);
    [t, x] = ode23s(KPC, tspan, x0, options);
    
    max_CN_values(i) = max(x(:,8));
    max_CN_valuesSS(i) = x(end,8);
    RespuestaSS(i) = (phi^5)*(x(end,3)+x(end,4)+x(end,5)+x(end,6)+x(end,7)+x(end,8));
    %max_x7_values(i) = x(end,7);
end

EmaxSim = max_CN_valuesSS(end);

figure;
semilogx(x0_values, max_CN_values, '-o');
xlabel('Total ligands');
ylabel('Maximal response');

figure;
semilogx(x0_values, max_CN_valuesSS, '-o');
xlabel('Total ligands');
ylabel('Maximal response');


figure;
semilogx(x0_values, max_CN_valuesSS, 'o', ...
    'MarkerSize', 4, ...                % Tamaño de marcadores un poco mayor
    'MarkerEdgeColor', 'k');
hold on;
semilogx(x0_values, RespuestaSS, '-', ...
    'LineWidth', 1, ...
    'Color', [0 0.4470 0.7410]);        % Mismo color azul para la línea
xlabel('$X_T$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('$\widehat R$ (response at steady-state)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Helvetica');
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
box off;
% --- Líneas horizontales E_max y E_max/2 ---
%[~, idx_max] = max(max_Tp_values);
%maxVal = max_Tp_values(idx_max);
maxVal = Emax;
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
x_half = interp1(max_CN_valuesSS, x0_values, half_val, 'spline');
hLineV = xline(x_half, ...
    'Color', [0.9290 0.6940 0.1250], ...  % naranja claro
    'LineStyle', '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '$EC_{50}$');
ylim([0 1.05 * maxVal]); 
legend([hLineH, hLineV, hLineH1], 'Interpreter', 'latex', 'Location', 'best');
hold off;


function dx = ODEKPRmcK(t, x, p)
    dx = zeros(4,1);
    dx(1) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4) + p(2) * x(5) + p(2) * x(6) + p(2) * x(7) + p(2) * x(8); %L
    dx(2) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4) + p(2) * x(5) + p(2) * x(6) + p(2) * x(7) + p(2) * x(8); %T
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3); %C0
    dx(4) = p(3) * x(3) - (p(2) + p(3)) * x(4); %C1  
    dx(5) = p(3) * x(4) - (p(2) + p(3)) * x(5); %C2  
    dx(6) = p(3) * x(5) - (p(2) + p(3)) * x(6); %C3  
    dx(7) = p(3) * x(6) - (p(2) + p(3)) * x(7); %C4  
    dx(8) = p(3) * x(7) - p(2) * x(8); %C5  
end