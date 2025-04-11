clear
clc

%%

load('resultadosCN_IFF.mat');

[~, idx_max] = max(CN_SS);
maxVal = CN_SS(idx_max);
half_val = maxVal / 2;

mitadInf = CN_SS(1:idx_max);
XT_inf   = XT_values(1:idx_max);

mitadSup = CN_SS(idx_max:end);
XT_sup   = XT_values(idx_max:end);

x_half = interp1(mitadInf, XT_inf, half_val, 'spline');
x_half2 = interp1(mitadSup, XT_sup, half_val, 'spline');

disp('Soluciones:');
disp([x_half, x_half2]);
