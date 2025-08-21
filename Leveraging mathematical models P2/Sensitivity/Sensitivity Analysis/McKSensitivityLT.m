clear

%% VALORES INICIALES INTEGRACIÓN DEL MODELO
    % initial values
    x0 = complex([100, 2e4, 0 , 0], 0); 
    % step size and time interval in days
    d = 1.0e-16; 
    tspan = 0.0:0.05:50;
    % kon koff
    p = complex([5e-5, 0.01, 1], 0);
    solution = sensitivity(x0, p, d, tspan); 

% 
%     % --------------- FORWARD RATE -----------------------------
% Vector de valores de koff
LTvect = logspace(0, log10(2e4), 2000);

h = waitbar(0,'Calculando...');
% Resultados con el número de filas de koff y en cada columna el instante
% temporal
results_matrix = zeros(length(LTvect), length(solution{1}(:, 1))); 
for i = 1:length(LTvect)

    x0 = complex([100, LTvect(i), 0 , 0], 0); 

    solution = sensitivity(x0, p, d, tspan);

    % COJO LA RESPUESTA QUE ME INTERESA:
    SolResponse = solution{4}(:, 3); 
    % Normalización de la respuesta
    newSol = (SolResponse .* LTvect(i)) ./ solution{4}(:, 1); 

    % En la fila que define un valor de koff
    results_matrix(i, :) = newSol;
    waitbar(i/length(LTvect), h);
end

close(h);

save('KPR1_LT.mat','tspan','LTvect','results_matrix');

results_matrix1 = log10(results_matrix);
LTvect1 = log10(LTvect);

figure('Position',[100 100 600 400]);
contourf(tspan, LTvect1, results_matrix, 10,'LineColor','k');
colormap(gray); colorbar
xlabel('Time (s)', 'FontSize',18);
ylabel('Total ligand (L_T)', 'FontSize',18);
title ('KPR-1', 'FontSize',18,'FontWeight','bold');
ticks_real = [1, 10, 100, 1e3, 1e4];
set(gca, 'YTick', log10(ticks_real));
set(gca, 'YTickLabel', {'1', '10', '100', '10^{3}', '10^{4}'});
set(gca, 'YDir', 'normal', 'FontSize', 16);

% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, LTvect, results_matrix);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Number of ligands', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('LTvect', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on



%% SOLUCION

% solution = sensitivity(x0, p, d, tspan); 
% 
% % solution{estado}(:, nºparametro)
% NewSolR = solution{4}(:, 1);
% 
% figure
% % Crear el gráfico
% plot(tspan, NewSolR);
% xlabel('t');
% legend;
% title('Sensitivity');
% 
% figure
% % COMPROBACIÓN
% neg = @(t,y)ODEKPRmcK(t, y, p);
% options = odeset('RelTol',1e-6,'AbsTol',1e-9);
% [t,x] = ode45(neg, tspan, x0, options);
% plot(t, x(:,4), 'DisplayName', 'x1');



%% FUNCIONES
function solution = sensitivity(x0, p, d, tspan)

    ST = @(t,y)ODEKPRmcK(t, y, p);
    options = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Refine', 1);
    [t,x] = ode23s(ST, tspan, x0, options);
    
    lp = length(p); ls = size(x, 1); lx = length(x0);
    % Crea un array de celdas de 1 fila y lx columnas. Cada celda puede contener datos de cualquier tipo, en este caso, matrices de ceros.
    solution = cell(1, lx);
    % Para cada índice i, se asigna una matriz de ceros de tamaño ls x (lp + 1) a la celda solution{i
    for i = 1:lx
        solution{i} = zeros(ls, lx + 1);
    end
    % El bucle itera sobre cada especie j y almacena la solución correspondiente en la primera columna de la 
    % matriz solution[j] dentro del diccionario solution. La primera columna se utiliza para almacenar la solución original (sin perturbaciones en los parámetros).
    for j = 1:lx
        solution{j}(:, 1) = x(:, j);
    end

    for j = 1:lx
        % La técnica de diferencias finitas complejas implica agregar una pequeña perturbación imaginaria a un parámetro 
        % para calcular la derivada parcial de la solución con respecto a ese parámetro.
        x0(j) = x0(j) + d * 1i; % Perturba el parámetro
        
        options = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Refine', 1);
        ST = @(t,y)ODEKPRmcK(t, y, p);
        [t,x] = ode23s(ST, tspan, x0, options);
        
        % Está destinada a restablecer el parámetro p[j] a su valor original, eliminando cualquier componente imaginaria que se haya agregado durante el proceso de perturbación.
        x0(j) = complex(real(x0(j)), 0);
        %  Toma la parte imaginaria de cada elemento en sol, donde sol es la matriz de soluciones del sistema ODE después de 
        % perturbar el parámetro correspondiente con una pequeña cantidad imaginaria d * im.
        % Divide la parte imaginaria de sol por d. Esto proporciona una aproximación de la derivada parcial de la solución 
        % con respecto al parámetro perturbado, utilizando diferencias finitas complejas.
        xSens = imag(x) ./ d;
        
        % Selecciona todas las filas y la columna j + 1 de la matriz solution[k]. 
        % La columna j + 1 se utiliza para almacenar la sensibilidad respecto al parámetro p[j] (perturbado).
        % Selecciona la fila k de la matriz sol (x), que contiene las sensibilidades calculadas para la especie k en todos los tiempos de evaluación.
        for k = 1:lx
            solution{k}(:, j + 1) = xSens(:, k);
        end


    end

    

end

function dx = ODEKPRmcK(t, x, p)
    dx = zeros(4,1);
    dx(1) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(2) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);
    dx(4) = p(3) * x(3) - p(2) * x(4);
end