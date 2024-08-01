clear

%% NEGATIVE II
% initial values
x0 = complex([0, 0, 0, 0, 0, 0, 0], 0); 
% step size and time interval in days
d = 1.0e-11; 
tspan = [0.0 20];
p = complex([1e-4, 10000, 3e4, 0.09, 10, 1.5, 0.04, 1.2e-6, 6e5, 1, 500, 10^4], 0);

%% McKeithan
x0 = complex([100, 2e4, 0, 0], 0);
tspan = [0.0 50];
d = 1.0e-16; 
p = complex([5e-5, 0.01, 1], 0);

%% SOLUCION

solution = sensitivity(x0, p, d, tspan); 

%% Plotear

% solution{estado}(:, nºparametro)

NewSolR = solution{4}(:, 1);

% Crear el gráfico
plot(NewSolR, 'DisplayName', 'x1');
xlabel('t');
ylabel('S');
legend;
title('Solución de la tercera especie sin perturbaciones');


%% FUNCIONES
function solution = sensitivity(x0, p, d, tspan)

    neg = @(t,y)ODEKPRmcK(t, y, p);
    options = odeset('RelTol',1e-3,'AbsTol',1e-3);
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
        
        options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        neg = @(t,y)ODEKPRmcK(t, y, p);
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
    dx = zeros(7,1);

    dx(1) = p(1) * (p(2) - x(1) - x(2) - x(3)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(5)) + p(4)) * x(1) + (p(7) + p(8)*x(4))*x(2);
    dx(2) = p(4)*x(1) + (p(7) + p(8)*x(4))*x(3) - ((1/p(5)) + p(4) + p(7) + p(8)*x(4))*x(2);
    dx(3) = p(4)*x(2) - ((1/p(5)) + p(7) + p(8)*x(4))*x(3);
    dx(4) = p(10)*(x(2) + x(6))*(p(9) - x(4)) - p(11)*x(4);
    dx(5) = p(1) * (p(12) - x(5) - x(6) - x(7)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(6)) + p(4)) * x(5) + (p(7) + p(8)*x(4))*x(6);
    dx(6) = p(4)*x(5) + (p(7) + p(8)*x(4))*x(7) - ((1/p(6)) + p(4) + p(7) + p(8)*x(4))*x(6);
    dx(7) = p(4)*x(6) - ((1/p(6)) + p(7) + p(8)*x(4))*x(7);

end

function dx = ODEKPRmcK(t, x, p)
    dx = zeros(4,1);
    dx(1) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(2) = - p(1) * x(1) * x(2) + p(2) * x(3) + p(2) * x(4);
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3);
    dx(4) = p(3) * x(3) - p(2) * x(4);
end