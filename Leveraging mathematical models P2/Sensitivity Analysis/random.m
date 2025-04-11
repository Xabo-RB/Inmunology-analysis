clear
% Parámetros
N = 5;       % Punto central donde ocurre el cambio
k = 0.5;     % Parámetro de crecimiento exponencial

% Valores de i
i_values = 0:0.1:5;  % Desde 0 hasta 10, con incrementos de 0.1

% Cálculo de alpha_i
alpha_values = exp(-k * (N - i_values));
%alpha_values = exp(k * (i_values));

% Graficar
figure;
plot(i_values, alpha_values, 'g', 'LineWidth', 1.5); % Curva en color verde
hold on;
xlabel('i', 'FontSize', 12);
ylabel('\alpha_i', 'FontSize', 12);
