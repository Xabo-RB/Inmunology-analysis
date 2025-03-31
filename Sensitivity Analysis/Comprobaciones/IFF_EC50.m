clear
clc

%%

load('resultadosCN_IFF.mat');  % Esto carga XT_values y CN_SS en el workspace

x_solutions = invert_interpolation(CN_SS, XT_values, half_val);
disp(x_solutions)

function x_vals = invert_interpolation(x_data, y_data, y_target)
    % Interpola x en función de y usando spline y encuentra todas las x para un y dado
    
    % Crear una interpolación densa
    y_dense = linspace(min(y_data), max(y_data), 1000);
    x_dense = interp1(y_data, x_data, y_dense, 'spline');
    
    % Calcular la diferencia respecto al valor deseado
    f = @(y) interp1(y_data, x_data, y, 'spline') - y_target;
    
    % Buscar cruces
    yy = y_dense;
    ff = x_dense - y_target;
    
    crossings = find(diff(sign(ff)) ~= 0);
    
    x_vals = [];
    for i = 1:length(crossings)
        y1 = yy(crossings(i));
        y2 = yy(crossings(i)+1);
        % Buscar el valor de y en que x = y_target
        y_root = fzero(f, [y1, y2]);
        % Interpolar x para ese y_root
        x_val = interp1(y_data, x_data, y_root, 'spline');
        x_vals(end+1) = x_val;
    end
end
