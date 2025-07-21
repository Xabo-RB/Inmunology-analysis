clear

%% VALORES INICIALES INTEGRACIÓN DEL MODELO
    %% Serial Triggering
    % initial values
    x0 = complex([0, 0, 0], 0); 
    % step size and time interval in days
    d = 1.0e-16; 
    tspan = 0.0:0.05:50;
    % p(1) = lambda, p(2) = phi, p(3) = s, p(4) = keff, p(5) = h, p(6) = L, p(7) = ki
    %p = complex([0.61, 0.0055, 0.0011, 0.001, 5, 1, 0.094], 0);

    % p(1) = varphi, p(2) = lambda, p(3) = sigma, p(4) = k_eff, p(5) = h, p(6) = k_i, p(7) = L
    p = complex([0.0055, 0.61,  0.0011, 0.001, 5, 0.094, 2e4], 0);

    solution = sensitivity(x0, p, d, tspan); 

%     % --------------- LT -----------------------------
LTvect = logspace(0, log10(2e4), 1000);

% Resultados con el número de filas de koff y en cada columna el instante
% temporal
results_matrix = zeros(length(LTvect), length(solution{3}(:, 1))); 
for i = 1:length(LTvect)

    p = complex([0.0055, 0.61,  0.0011, 0.001, 5, 0.094, LTvect(i)], 0);

    solution = sensitivity(x0, p, d, tspan);

    % COJO LA RESPUESTA QUE ME INTERESA:
    SolResponse = solution{3}(:, 8); 
    % Normalización de la respuesta
    newSol = (SolResponse .* LTvect(i)) ./ solution{3}(:, 1); 

    % En la fila que define un valor de koff
    results_matrix(i, :) = newSol;
end


results_matrix1 = log10(abs(results_matrix));
LTvect1 = log10(LTvect);

figure('Position',[100 100 600 400]);
contourf(tspan, LTvect1, results_matrix, 10,'LineColor','k');
colormap(gray); colorbar
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel('Total ligand (L_T)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
title ('ST', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
ticks_real = [1, 10, 100, 1e3, 1e4];
set(gca, 'YTick', log10(ticks_real));
set(gca, 'YTickLabel', {'1', '10', '100', '10^{3}', '10^{4}'});


%% FUNCIONES
function solution = sensitivity(x0, p, d, tspan)

    ST = @(t,y)ODESTMod(t, y, p);
    options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
    [t,x] = ode23s(ST, tspan, x0, options);
    
    lp = length(p); ls = size(x, 1); lx = length(x0);
    % Crea un array de celdas de 1 fila y lx columnas. Cada celda puede contener datos de cualquier tipo, en este caso, matrices de ceros.
    solution = cell(1, lx);
    % Para cada índice i, se asigna una matriz de ceros de tamaño ls x (lp + 1) a la celda solution{i
    for i = 1:lx
        solution{i} = zeros(ls, lp + 1);
    end
    % El bucle itera sobre cada especie j y almacena la solución correspondiente en la primera columna de la 
    % matriz solution[j] dentro del diccionario solution. La primera columna se utiliza para almacenar la solución original (sin perturbaciones en los parámetros).
    for j = 1:lx
        solution{j}(:, 1) = x(:, j);
    end

    for j = 1:lp
        % La técnica de diferencias finitas complejas implica agregar una pequeña perturbación imaginaria a un parámetro 
        % para calcular la derivada parcial de la solución con respecto a ese parámetro.
        p(j) = p(j) + d * 1i; % Perturba el parámetro
        
        options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
        ST = @(t,y)ODESTMod(t, y, p);
        [t,x] = ode23s(ST, tspan, x0, options);
        
        % Está destinada a restablecer el parámetro p[j] a su valor original, eliminando cualquier componente imaginaria que se haya agregado durante el proceso de perturbación.
        p(j) = complex(real(p(j)), 0);
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

function dx = ODESerialTrig(t, x, p)
    dx = zeros(4,1);
    dx(1) = -p(1) * p(2) * (x(1) - x(2)) + p(3) * (1 - x(1));
    dx(2) = p(2) * (x(1) - x(2)) + p(3) * (1 - x(2)) - p(4) * (x(2)^p(5)) * (p(6)^p(5));
    dx(3) = p(4) * (x(2)^p(5)) * (p(6)^p(5)) - p(7) * x(3);
    dx(4) = x(1) / (p(1) + 1) + ((x(2) + x(3)) * p(1) / (p(1) + 1));
end

function dx = ODESTMod(t, x, p)
    dx = zeros(3, 1);

    % Definir las ecuaciones diferenciales
    % p(1) = varphi, p(2) = lambda, p(3) = sigma, p(4) = k_eff, p(5) = h, p(6) = k_i, p(7) = L
    dx(1) = -p(1) * (x(1) - (1 / p(2)) * x(2)) + p(3) * (1 / (1 + p(2)) - x(1));  % Ecuación para dS/dt
    dx(2) = p(1) * (x(1) - (1 / p(2)) * x(2)) + p(3) * (p(2) / (1 + p(2)) - x(2)) - p(4) * x(2)^p(5) * p(7)^p(5);  % Ecuación para dT/dt
    dx(3) = p(4) * x(2)^p(5) * p(7)^p(5) - p(6) * x(3);  % Ecuación para dA/dt
end


