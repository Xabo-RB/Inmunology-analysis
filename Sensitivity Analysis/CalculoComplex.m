clear

%% VALORES INICIALES INTEGRACIÓN DEL MODELO
    %% NEGATIVE II
    % initial values
    x0 = complex([0, 0, 0, 0, 0, 0, 0, 0], 0); 
    % step size and time interval in days
    d = 1.0e-16; 
    tspan = 0.0:0.05:20;
    %k = p[1] = kon,  L1 = p[2], R = p[3], phi = p[4], tao1 = p[5], 
    % tao2 = p[6], b = p[7], gamma = p[8], ST = p[9], alpha = p[10], 
    % beta = p[11], L2 = p[12]
    p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, 0.04, 1.2e-6, 6e5, 1, 500, 10^4], 0);
    solution = sensitivity(x0, p, d, tspan); 

    % --------------- TAO1 / KOFF1 -----------------------------
% Vector de valores de koff
koffVect = 0.001:0.001:1;
tao1Vect = 1./koffVect;
% Resultados con el número de filas de koff y en cada columna el instante
% temporal
results_matrix = zeros(length(tao1Vect), length(solution{8}(:, 1))); 
for i = 1:length(tao1Vect)

    p = complex([1e-4, 10000, 3e4, 0.09, tao1Vect(i), 1.5, 0.04, 1.2e-6, 6e5, 1, 500, 10^4], 0);

    solution = sensitivity(x0, p, d, tspan);

    % COJO LA RESPUESTA QUE ME INTERESA:
    SolResponse = solution{8}(:, 6); 
    % Normalización de la respuesta
    newSol = (SolResponse .* tao1Vect(i)) ./ solution{8}(:, 1); 

    % En la fila que define un valor de koff
    results_matrix(i, :) = newSol;
end

% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, tao1Vect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Dissociate rate');
% title('Negative regulator II');
% set(gca, 'YDir', 'normal');

inferno = csvread('inferno_colormap.csv');
figure('Position', [100, 100, 600, 400]);
imagesc(tspan, koffVect, results_matrix);
colormap(inferno);
cb = colorbar;
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel('Unbinding rate', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YDir', 'normal');
hold on

    % --------------- TAO1 / KOFF1 para S(t) -----------------------------
% % Vector de valores de koff
% koffVect = 0.001:0.001:1;
% tao1Vect = 1./koffVect;
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(tao1Vect), length(solution{4}(:, 1))); 
% for i = 1:length(tao1Vect)
% 
%     p = complex([1e-4, 10000, 3e4, 0.09, tao1Vect(i), 1.5, 0.04, 1.2e-6, 6e5, 1, 500, 10^4], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{4}(:, 6); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* tao1Vect(i)) ./ solution{4}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, tao1Vect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Dissociate rate');
% title('Quantity of active SHP-1 molecules');
% set(gca, 'YDir', 'normal');

    % --------------- b (dephosphorylation rate) SHP1 -----------------------------
% % Vector de valores de koff
% bVect = 0.004:0.001:0.44;
% bVect = linspace(0.004, 0.44, 1000); % 1000% por arriba y por abajo
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(bVect), length(solution{4}(:, 1))); 
% for i = 1:length(bVect)
% 
%     p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, bVect(i), 1.2e-6, 6e5, 1, 500, 10^4], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{4}(:, 8); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* bVect(i)) ./ solution{4}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, bVect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Spontaneous dephosphorylation rate');
% title('Sensitivity of [SHP-1] to $b$', 'Interpreter', 'latex');
% set(gca, 'YDir', 'normal');
% hold on

    % --------------- gamma SHP1-----------------------------
% % Vector de valores de koff
% gammaVect = 1e-8:1e-7:1e-4;
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(gammaVect), length(solution{4}(:, 1))); 
% for i = 1:length(gammaVect)
% 
%     p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, 0.04, gammaVect(i), 6e5, 1, 500, 10^4], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{4}(:, 9); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* gammaVect(i)) ./ solution{4}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, gammaVect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Dephosphorylation rate by SHP-1');
% title('Sensitivity of [SHP-1] to $\gamma$', 'Interpreter', 'latex');
% set(gca, 'YDir', 'normal');
% hold on


    % --------------- b (dephosphorylation rate) -----------------------------
% % Vector de valores de koff
% bVect = 0.004:0.001:0.44;
% bVect = linspace(0.004, 0.44, 1000); % 1000% por arriba y por abajo
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(bVect), length(solution{8}(:, 1))); 
% for i = 1:length(bVect)
% 
%     p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, bVect(i), 1.2e-6, 6e5, 1, 500, 10^4], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{8}(:, 8); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* bVect(i)) ./ solution{8}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, bVect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Spontaneous dephosphorylation rate');
% title('Sensitivity to $b$', 'Interpreter', 'latex');
% set(gca, 'YDir', 'normal');
% hold on

    % --------------- gamma -----------------------------
% % Vector de valores de koff
% gammaVect = 1e-8:1e-7:1e-4;
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(gammaVect), length(solution{8}(:, 1))); 
% for i = 1:length(gammaVect)
% 
%     p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, 0.04, gammaVect(i), 6e5, 1, 500, 10^4], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{8}(:, 9); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* gammaVect(i)) ./ solution{8}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% %inferno = flipud(inferno);
% figure; 
% %imagesc(tspan, koffVect, results_matrix); 
% imagesc(tspan, gammaVect, results_matrix); 
% colormap(inferno);
% cb = colorbar;
% cb.Label.String = 'Sensitivity';
% xlabel('Time (s)');
% ylabel('Dephosphorylation rate by SHP-1');
% title('Sensitivity to $\gamma$', 'Interpreter', 'latex');
% set(gca, 'YDir', 'normal');
% hold on

    %% McKeithan
%     x0 = complex([100, 2e4, 0, 0], 0);
%     % Fijar el paso del intervalo de tiempo para que no varíe el tamaño del
%     % vector con cada resultado
%     tspan = 0:0.1:50;
%     d = 1.0e-16; 
%     p = complex([5e-5, 0.01, 1], 0);
% 
% % Vector de valores de koff
% koffVect = 0.001:0.001:1;
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(koffVect), length(solution{4}(:, 1))); 
% for i = 1:length(koffVect)
%     p = complex([5e-5, koffVect(i), 1], 0);
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{4}(:, 3); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* koffVect(i)) ./ solution{4}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% inferno = csvread('inferno_colormap.csv');
% figure; 
% imagesc(tspan, koffVect, results_matrix); 
% colormap(inferno);
% colorbar;
% xlabel('Time (s)');
% ylabel('Dissociate rate (koff)');
% title('McKeithan');
% set(gca, 'YDir', 'normal');


%% SOLUCION

solution = sensitivity(x0, p, d, tspan); 

% solution{estado}(:, nºparametro)
NewSolR = solution{8}(:, 1);

figure
% Crear el gráfico
plot(tspan, NewSolR);
xlabel('t');
legend;
title('Sensitivity');

% % COMPROBACIÓN
% neg = @(t,y)ODEKPRmcK(t, y, p);
% options = odeset('RelTol',1e-6,'AbsTol',1e-9);
% [t,x] = ode45(neg, tspan, x0, options);
% plot(t, x(:,4), 'DisplayName', 'x1');



%% FUNCIONES
function solution = sensitivity(x0, p, d, tspan)

    neg = @(t,y)negativeII(t, y, p);
    options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
    [t,x] = ode45(neg, tspan, x0, options);
    
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
        neg = @(t,y)negativeII(t, y, p);
        [t,x] = ode45(neg, tspan, x0, options);
        
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

function dx = negativeII(t, x, p)
    dx = zeros(8,1);

    dx(1) = p(1) * (p(2) - x(1) - x(2) - x(3)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(5)) + p(4)) * x(1) + (p(7) + p(8)*x(4))*x(2);
    dx(2) = p(4)*x(1) + (p(7) + p(8)*x(4))*x(3) - ((1/p(5)) + p(4) + p(7) + p(8)*x(4))*x(2);
    dx(3) = p(4)*x(2) - ((1/p(5)) + p(7) + p(8)*x(4))*x(3);
    dx(4) = p(10)*(x(2) + x(6))*(p(9) - x(4)) - p(11)*x(4);
    dx(5) = p(1) * (p(12) - x(5) - x(6) - x(7)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(6)) + p(4)) * x(5) + (p(7) + p(8)*x(4))*x(6);
    dx(6) = p(4)*x(5) + (p(7) + p(8)*x(4))*x(7) - ((1/p(6)) + p(4) + p(7) + p(8)*x(4))*x(6);
    dx(7) = p(4)*x(6) - ((1/p(6)) + p(7) + p(8)*x(4))*x(7);
    dx(8) = p(4)*x(2) - ((1/p(5)) + p(7) + p(8)*x(4))*x(3) + p(4)*x(6) - ((1/p(6)) + p(7) + p(8)*x(4))*x(7);

end

function dx = ODEKPRmcK(t, x, p)
    dx = zeros(4,1);
    dx(1) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(2) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);
    dx(4) = p(3) * x(3) - p(2) * x(4);
end