clear

%% VALORES INICIALES INTEGRACIÓN DEL MODELO
    % initial values
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0, 0], 0); 
    % step size and time interval in days
    d = 1.0e-16; 
    tspan = 0.0:0.05:200;
    % kon koff
    p = complex([5e-5, 0.01, 1, 4.4e-4, 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
    solution = sensitivity(x0, p, d, tspan); 

    % --------------- UNBINDING RATE -----------------------------
% Vector de valores de koff
% koffVect = 0.001:0.001:1;
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(koffVect), length(solution{1}(:, 1))); 
% for i = 1:length(koffVect)
% 
%     p = complex([5e-5, koffVect(i), 1, 4.4e-4, 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{10}(:, 3); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* koffVect(i)) ./ solution{10}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% save('nf2_KOFF.mat','tspan','koffVect','results_matrix');
% 
% results_matrix1 = abs(results_matrix);
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, koffVect, results_matrix, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% % results_matrix1 = abs(log10(results_matrix));
% % 
% % figure('Position', [100, 100, 600, 400]);
% % contourf(tspan, koffVect, results_matrix1, 10, 'LineColor', 'k');
% % colormap(gray);
% % colorbar;
% % xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% % title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% % set(gca, 'YDir', 'normal');
% % hold on
% 
% results_matrix2 = log10(abs(results_matrix));
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, koffVect, results_matrix, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% %ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% ylabel('k_{off}', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'tex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
% hold on
% 
% % 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, koffVect, results_matrix2, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, koffVect, results_matrix2);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% % hold on


%     % --------------- BINDING RATE -----------------------------
% Vector de valores de koff
% konVect = 4e-6:1e-6:2e-2;
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% h = waitbar(0,'Calculando...');
% results_matrix = zeros(length(konVect), length(solution{1}(:, 1))); 
% for i = 1:length(konVect)
% 
%     p = complex([konVect(i), 0.01, 1, 4.4e-4, 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{10}(:, 2); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* konVect(i)) ./ solution{10}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
%     waitbar(i/length(konVect), h);
% end
% 
% close(h);
% 
% save('nf2_Kon.mat','tspan','konVect','results_matrix');
% 
% results_matrix1 = log10(abs(results_matrix));
% 
% fig = figure('Visible','off');
% [C,h] = contourf(tspan, konVect, results_matrix1, 10,'LineColor','k');
% levels = h.LevelList;
% fprintf('levels = [%s]\n', num2str(levels, '%.4g, '));
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, konVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('k_{on}', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'tex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
% ylim([min(konVect) max(konVect)]);
% yticks = 10.^(ceil(log10(min(konVect))):floor(log10(max(konVect))));
% set(gca, 'YTick', yticks);
% set(gca, 'YTickLabel', compose('10^{%d}', log10(yticks)));
% hold on
% 
% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, konVect, results_matrix);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Binding rate', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on

% 
%     % --------------- FORWARD RATE -----------------------------
%Vector de valores de koff
% kpVect = linspace(0.001, 10, 2000); 
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% h = waitbar(0,'Calculando...');
% results_matrix = zeros(length(kpVect), length(solution{1}(:, 1))); 
% for i = 1:length(kpVect)
% 
%     p = complex([5e-5, 0.01, kpVect(i), 4.4e-4, 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{10}(:, 4); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* kpVect(i)) ./ solution{10}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
%     waitbar(i/length(kpVect), h);
% end
% 
% close(h);
% 
% save('KPR_NF2_kp.mat','tspan','kpVect','results_matrix');
% 
% results_matrix1 = log10(abs(results_matrix));
% 
% fig = figure('Visible','off');
% [C,h] = contourf(tspan, kpVect, results_matrix1, 10,'LineColor','k');
% levels = h.LevelList;
% fprintf('levels = [%s]\n', num2str(levels, '%.4g, '));
% 
% results_matrix1 = log10(abs(results_matrix));
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, kpVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('k_{p}', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'tex');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
% hold on
% 
% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, kpVect, results_matrix);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Phosphorylation rate', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on

% 
% %    --------------- b (dephosphorylation rate) TO r(T) -----------------------------
% % Vector de valores de koff
% bVect = 0.004:0.001:0.44;
% bVect = linspace(0.05, 0.1, 500); % 1000% por arriba y por abajo
% bVect = linspace(1e-3, 1, 1000); 
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(bVect), length(solution{1}(:, 1))); 
% h = waitbar(0,'Calculando...');
% for i = 1:length(bVect)
% 
%     p = complex([5e-5, 0.01, 1, 4.4e-4, bVect(i), 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{10}(:, 6); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* bVect(i)) ./ solution{10}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
%     waitbar(i/length(bVect), h);
% end
% 
% close(h);
% 
% results_matrix1 = log10(abs(results_matrix));
% 
% save('NF2_Rto_b.mat','tspan','bVect','results_matrix');
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, bVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Spontaneous dephospho. ($b$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of $R$ to parameter $b$ (2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% % inferno = csvread('inferno_colormap.csv');
% % figure('Position', [100, 100, 600, 400]);
% % imagesc(tspan, bVect, results_matrix);
% % colormap(inferno);
% % cb = colorbar;
% % xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % ylabel('Phosphorylation rate', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % title('Sensitivity of $R$ to parameter $b$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% % set(gca, 'YDir', 'normal');
% % hold on
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, bVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Spontaneous dephospho. ($b$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of $R$ to $b$ (2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, bVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel({'Parameter b','Spontaneous dephospho.'}, ...
%        'Interpreter','tex', 'FontSize',18, 'Color','k', 'FontWeight','normal');
% title('Sensitivity of R to b', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
% hold on
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, bVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel({'Parameter b','Spontaneous dephospho.'}, ...
%        'Interpreter','tex', 'FontSize',18, 'Color','k', 'FontWeight','normal');
% title('Sensitivity of R to b', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
% hold on

% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, bVect, results_matrix);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Spontaneous dephospho. ($b$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of [SHP-1] to $b$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on


    % --------------- gamma R-----------------------------
% Vector de valores de koff
gammaVect = 1e-8:1e-7:1e-4;
gammaVect = linspace(1e-5,1e-1,1000);

% Resultados con el número de filas de koff y en cada columna el instante
% temporal
results_matrix = zeros(length(gammaVect), length(solution{1}(:, 1))); 
h = waitbar(0,'Calculando...');
for i = 1:length(gammaVect)

    p = complex([5e-5, 0.01, 1, gammaVect(i), 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);

    solution = sensitivity(x0, p, d, tspan);

    % COJO LA RESPUESTA QUE ME INTERESA:
    SolResponse = solution{10}(:, 5); 
    % Normalización de la respuesta
    newSol = (SolResponse .* gammaVect(i)) ./ solution{10}(:, 1); 

    % En la fila que define un valor de koff
    results_matrix(i, :) = newSol;
    waitbar(i/length(gammaVect), h);
end

close(h);

save('NF2_Rto_gamma.mat','tspan','gammaVect','results_matrix');

results_matrix1 = log10(abs(results_matrix));

figure('Position', [100, 100, 600, 400]);
contourf(tspan, gammaVect, results_matrix, 10, 'LineColor', 'k');
colormap(gray);
colorbar;
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
title('Sensitivity of $R$ to parameter $\gamma$ (KPR-NF2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YDir', 'normal');
hold on

figure('Position', [100, 100, 600, 400]);
contourf(tspan, gammaVect, results_matrix1, 10, 'LineColor', 'k');
colormap(gray);
colorbar;
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
title('Sensitivity of $R$ to parameter $\gamma$ (KPR-NF2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YDir', 'normal');
hold on

inferno = csvread('inferno_colormap.csv');
figure('Position', [100, 100, 600, 400]);
imagesc(tspan, gammaVect, results_matrix);
colormap(inferno);
cb = colorbar;
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
title('Sensitivity of $R$ to parameter $\gamma$ (KPR-NF2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YDir', 'normal');
hold on

figure('Position', [100, 100, 600, 400]);
contourf(tspan, gammaVect, results_matrix1, 10, 'LineColor', 'k');
colormap(gray);
colorbar;
xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
ylabel({'Parameter \gamma','Dephospho. by SHP-1'}, ...
       'Interpreter','tex', 'FontSize',18, 'Color','k', 'FontWeight','normal');
title('Sensitivity of R to \gamma', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YDir', 'normal', 'FontSize', 16, 'YScale', 'log');
hold on


    % --------------- gamma SHP1-----------------------------
% % Vector de valores de koff
% gammaVect = 1e-8:1e-7:1e-4;
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(gammaVect), length(solution{1}(:, 1))); 
% for i = 1:length(gammaVect)
% 
%     p = complex([5e-5, 0.01, 1, gammaVect(i), 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{9}(:, 5); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* gammaVect(i)) ./ solution{9}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% save('NF2_Rto_gammatoS.mat','tspan','gammaVect','results_matrix');
% 
% results_matrix1 = log10(abs(results_matrix));
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, gammaVect, results_matrix, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of [SHP-1] to $\gamma$ (KPR-NF2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, gammaVect, results_matrix1, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of [SHP-1] to $\gamma$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on

% inferno = csvread('inferno_colormap.csv');
% figure('Position', [100, 100, 600, 400]);
% imagesc(tspan, gammaVect, results_matrix);
% colormap(inferno);
% cb = colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('Dephospho. rate by SHP-1 ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% title('Sensitivity of [SHP-1] to $\gamma$ (KPR-NF2)', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on

    % --------------- UNBINDING RATE -----------------------------
% % Vector de valores de koff
% koffVect = 0.001:0.001:1;
% 
% % Resultados con el número de filas de koff y en cada columna el instante
% % temporal
% results_matrix = zeros(length(koffVect), length(solution{1}(:, 1))); 
% for i = 1:length(koffVect)
% 
%     p = complex([5e-5, koffVect(i), 1, 4.4e-4, 0.04, 1, 2e-4, 6e5, exp(-1*2), exp(-1*1), exp(0)], 0);
% 
%     solution = sensitivity(x0, p, d, tspan);
% 
%     % COJO LA RESPUESTA QUE ME INTERESA:
%     SolResponse = solution{9}(:, 3); 
%     % Normalización de la respuesta
%     newSol = (SolResponse .* koffVect(i)) ./ solution{9}(:, 1); 
% 
%     % En la fila que define un valor de koff
%     results_matrix(i, :) = newSol;
% end
% 
% save('NF2_kofftoS.mat','tspan','koffVect','results_matrix');
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, koffVect, results_matrix, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('Sensitivity of [SHP-1] to $k_{off}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% % results_matrix1 = abs(log10(results_matrix));
% % 
% % figure('Position', [100, 100, 600, 400]);
% % contourf(tspan, koffVect, results_matrix1, 10, 'LineColor', 'k');
% % colormap(gray);
% % colorbar;
% % xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % ylabel('Unbinding rate', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % title('KPR-NF2', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% % set(gca, 'YDir', 'normal');
% % hold on
% 
% results_matrix2 = log10(abs(results_matrix));
% 
% figure('Position', [100, 100, 600, 400]);
% contourf(tspan, koffVect, results_matrix2, 10, 'LineColor', 'k');
% colormap(gray);
% colorbar;
% xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('Sensitivity of [SHP-1] to $k_{off}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% set(gca, 'YDir', 'normal');
% hold on
% 
% % inferno = csvread('inferno_colormap.csv');
% % figure('Position', [100, 100, 600, 400]);
% % imagesc(tspan, koffVect, results_matrix2);
% % colormap(inferno);
% % cb = colorbar;
% % xlabel('Time (s)', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal');
% % ylabel('$k_{off}$', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'normal', 'Interpreter', 'latex');
% % title('Sensitivity of [SHP-1] to $k_{off}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
% % set(gca, 'YDir', 'normal');
% % hold on




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

% % COMPROBACIÓN
% neg = @(t,y)ODEKPRmcK(t, y, p);
% options = odeset('RelTol',1e-6,'AbsTol',1e-9);
% [t,x] = ode45(neg, tspan, x0, options);
% plot(t, x(:,4), 'DisplayName', 'x1');



%% FUNCIONES
function solution = sensitivity(x0, p, d, tspan)

    ST = @(t,y)ODEKPRNegFeed2(t, y, p);
    options = odeset('RelTol',1e-5,'AbsTol',1e-5, 'Refine', 1);
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
        
        options = odeset('RelTol',1e-5,'AbsTol',1e-5, 'Refine', 1);
        ST = @(t,y)ODEKPRNegFeed2(t, y, p);
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

function dx = ODEKPRNegFeed2(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(10, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8)); % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8)); % P
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(9)) * x(4); % C0
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(4) + (p(5) + p(4) * x(9)) * x(5); % C1
    dx(5) = p(3) * x(4) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(5) + (p(5) + p(4) * x(9)) * x(6); % C2
    dx(6) = p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7); % C3
    dx(7) = p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8); % C4
    dx(8) = p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8); % C5
    dx(9) = p(7) * x(4) * (p(8) - x(9)) - p(6) * x(9); % S
    dx(10) = p(9)*(p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7)) + ...
        p(10)*(p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8)) + ...
        p(11)*(p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8));
end
