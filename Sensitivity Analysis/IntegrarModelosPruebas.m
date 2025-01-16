clear

% %% Parameters
% kon =       5e-5;
% koff =      0.01;
% kp =        0.01;
% 
% %% Modelo con solvers de Matlab
% x0 = [100, 2e4, 0, 0];
% keithan = @(t,y)KPRMcK1(t,y,kon, koff, kp);
% options = odeset('RelTol',1e-3,'AbsTol',1e-3);
% [t,x] = ode45(keithan, [0 500], x0, options);
% 
% %% Graficar Resultados
% % figure;
% % plot(t, x);
% % xlabel('Tiempo');
% % ylabel('Concentraciones');
% % legend('P(t)', 'T(t)', 'C0(t)', 'C1(t)');
% % title('Solución del sistema de ODEs');
% 
% figure
% plot(t, x(:, 4));
% xlabel('Tiempo');
% ylabel('Concentraciones');
% legend('C1(t)');
% title('Solución del sistema de ODEs');

%% NEGATIVE II
%% Parameters
% k = 1e-4; phi = 0.09; 
% b = 0.04; gamma = 1.2e-6;
% tao1 = 10; tao2 = 1.5;
% alpha = 1; beta = 500; 
% L1 = 10000; R = 3e4; L2 = 10^4; ST = 6e5;
% %% Modelo con solvers de Matlab
% x0 = [0, 0, 0, 0, 0, 0, 0];
% negII = @(t,y)NegativeFeedbackII(t,y, k, L1, R, phi, tao1, tao2, b, gamma, ST, alpha, beta, L2);
% options = odeset('RelTol',1e-3,'AbsTol',1e-3);
% [t,x] = ode45(negII, [0 20], x0, options);
% 
% figure
% plot(t, x(:, 3));
% xlabel('Tiempo');
% ylabel('');
% legend('C1(t)');
% title('');

%% NEGATIVE I
%% Parameters
% Parámetros del modelo
kon = 5e-5;       % Constante de asociación
koff = 0.01;      % Constante de disociación
kp = 1;           % Constante de producción
gamma = 4.4e-4;   % Factor gamma
b = 0.04;         % Factor b
beta = 1;         % Factor beta
alpha = 2e-4;     % Factor alpha
ST = 6e5;         % Concentración máxima de S

%% Modelo con solvers de Matlab
x0 = [100, 2e4, 0, 0, 0, 0];
negII = @(t,y)NegativeFeedbackI(t, x, kon, koff, kp, b, gamma, ST, alpha, beta);
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,x] = ode45(negII, [0 200], x0, options);

figure
plot(t, x(:, 3));
xlabel('Tiempo');
ylabel('');
legend('C1(t)');
title('');



%% SERIAL TRIGGERING
%% Parameters
% lambda = 0.61; phi = 0.0055; s = 0.0011; ki = 0.094;
% keff = 0.001; 
% L = 1; h = 5;
% 
% %% Modelo con solvers de Matlab
% x0 = [0, 0, 0, 1];
% ST = @(t,y)SerialTriggering(t, y, lambda, phi, s, keff, ki, L, h);
% options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
% [t,x] = ode45(ST, [0 100], x0, options);
% 
% figure
% plot(t, x(:, 4));
% xlabel('Tiempo');
% ylabel('Concentraciones');
% legend('C1(t)');
% title('Solución del sistema de ODEs');

function dy = KPRMcK1(t, y, kon, koff, kp)
    
    dy = zeros(4,1);
    
    dy(1) = - kon * y(1) * y(2) + koff*y(3) + koff*y(4);       % P(t)
    dy(2) = - kon * y(1) * y(2) + koff*y(3) + koff*y(4);      % T(t)
    dy(3) = kon * y(1) * y(2) - (koff + kp)*y(3);              % C0(t)
    dy(4) = kp*y(3) - koff*y(4);                            % C1(t)
        
end

function dx = NegativeFeedbackII(t, x, k, L1, R, phi, tao1, tao2, b, gamma, ST, alpha, beta, L2)

    dx = zeros(7, 1);
    
    dx(1) = k * (L1 - x(1) - x(2) - x(3)) * (R - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/tao1) + phi) * x(1) + (b + gamma * x(4)) * x(2);
    dx(2) = phi * x(1) + (b + gamma * x(4)) * x(3) - ((1/tao1) + phi + b + gamma * x(4)) * x(2);
    dx(3) = phi * x(2) - ((1/tao1) + b + gamma * x(4)) * x(3);
    dx(4) = alpha * (x(2) + x(6)) * (ST - x(4)) - beta * x(4);
    dx(5) = k * (L2 - x(5) - x(6) - x(7)) * (R - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/tao2) + phi) * x(5) + (b + gamma * x(4)) * x(6);
    dx(6) = phi * x(5) + (b + gamma * x(4)) * x(7) - ((1/tao2) + phi + b + gamma * x(4)) * x(6);
    dx(7) = phi * x(6) - ((1/tao2) + b + gamma * x(4)) * x(7);

end

function dx = NegativeFeedbackI(t, x, kon, koff, kp, b, gamma, ST, alpha, beta)

    dx = zeros(6, 1);
    
    dx(1) = - kon * x(1) * x(2) + koff*x(3) + koff*x(4) + koff*x(5);
    dx(2) = - kon * x(1) * x(2) + koff*x(3) + koff*x(4) + koff*x(5);
    dx(3) = kon * x(1) * x(2) - (koff + kp)*x(3) + (b + gamma*x(6))*x(4);
    dx(4) = kp*x(3) - (koff + kp + b + gamma*x(6))*x(4) + (b + gamma*x(6))*x(5);
    dx(5) = kp*x(4) - (koff + b + gamma*x(6))*x(5);
    dx(6) = alpha*x(4)*(ST - x(6)) - beta*x(6);

end


function dx = SerialTriggering(t, x, lambda, phi, s, keff, ki, L, h)

    dx = zeros(4, 1);
    
    % S(t) = x(1)
    % T(t) = x(2)
    % A(t) = x(3)
    % Y(t) = x(4)

    dx(1) = -lambda * phi * (x(1) - x(2)) + s * (1 - x(1));
    dx(2) =  phi * (x(1) - x(2)) + s * (1 - x(2)) - keff * (x(2)^h) * (L^h);
    dx(3) =  keff * (x(2)^h) * (L^h) - ki * x(3);
    dx(4) =  x(1) / (lambda + 1) + ((x(2) + x(3)) * lambda / (lambda + 1));

end
