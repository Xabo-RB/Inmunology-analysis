clear

%% Parameters
kon =       5e-5;
koff =      0.01;
kp =        0.01;

%% Modelo con solvers de Matlab
x0 = [100, 2e4, 0, 0];
keithan = @(t,y)KPRMcK1(t,y,kon, koff, kp);
options = odeset('RelTol',1e-3,'AbsTol',1e-3);
[t,x] = ode45(keithan, [0 200], x0, options);

%% Graficar Resultados
figure;
plot(t, x);
xlabel('Tiempo');
ylabel('Concentraciones');
legend('P(t)', 'T(t)', 'C0(t)', 'C1(t)');
title('Solución del sistema de ODEs');

figure
plot(t, x(:, 4));
xlabel('Tiempo');
ylabel('Concentraciones');
legend('C1(t)');
title('Solución del sistema de ODEs');


function dy = KPRMcK1(t, y, kon, koff, kp)
    
    dy = zeros(4,1);
    
    dy(1) = - kon * y(1) * y(2) + koff*y(3) + koff*y(4);       % P(t)
    dy(2) = - kon * y(1) * y(2) + koff*y(3) + koff*y(4);      % T(t)
    dy(3) = kon * y(1) * y(2) - (koff + kp)*y(3);              % C0(t)
    dy(4) = kp*y(3) - koff*y(4);                            % C1(t)
        
end
