clear

% Parámetros
kon = 5e-5; % Valor ejemplo
koff = 0.01; % Valor ejemplo
kp = 1; % Valor ejemplo

P0 = 2e4; T0 = 100; C0_0 = 0; C1_0 = 0; C2_0 = 0; C3_0 = 0; C4_0 = 0;

y0 = [P0; T0; C0_0; C1_0; C2_0; C3_0; C4_0];

tspan = [0 100];

[t, y] = ode45(@(t, y) mck(t, y, kon, koff, kp), tspan, y0);

plot(t, y(:, 3:end));
%legend('P(t)', 'T(t)', 'C0(t)', 'C1(t)', 'C2(t)', 'C3(t)', 'C4(t)');
legend('C0(t)', 'C1(t)', 'C2(t)', 'C3(t)', 'C4(t)');
xlabel('Tiempo (t)');
ylabel('Unidades');
title('');

E_max = ((kp / (kp + koff))^4) * (y(:, 2) + y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6) + y(:, 7));
figure
plot(t,E_max)

TT = y(:, 2) + y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6) + y(:, 7);

clear

% Parámetros
k1 = 10;         
k3 = 1;       
kmenos1 = 0.1;  
k2 = 1;       
kmenos2 = 10;  
w = 1;        

% Valores iniciales
T0      = 100;        
X0      = 2;        
Tp0     = 0;        
D0      = 0;          
C0_0    = 0;     
C1_0    = 0;        
C2_0    = 0;        

% Vector de condiciones iniciales
y0 = [T0; X0; Tp0; D0; C0_0; C1_0; C2_0];

% Tiempo de simulación
tspan = [0 200]; % Tiempo inicial y final

% Resolver el sistema de ecuaciones diferenciales
[t, y] = ode45(@(t, y) kpc(t, y, k1, k3, kmenos1, k2, kmenos2, w), tspan, y0);

% Graficar resultados
plot(t, y);
legend('T(t)', 'X(t)', 'Tp(t)', 'D(t)', 'C0(t)', 'C1(t)', 'C2(t)');
xlabel('Tiempo (t)');
ylabel('Concentraciones');

TT = y(:, 1) + y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6) + y(:, 7);
XT = y(:, 2) + y(:, 4) + y(:, 5) + y(:, 6) + y(:, 7);
CT = y(:, 5) + y(:, 6) + y(:, 7);

% Función que define el sistema de ecuaciones diferenciales
function dydt = kpc(t, y, k1, k3, kmenos1, k2, kmenos2, w)
    % Variables
    T = y(1);
    X = y(2);
    Tp = y(3);
    D = y(4);
    C0 = y(5);
    C1 = y(6);
    C2 = y(7);

    % Ecuaciones diferenciales
    dTdt = -k1 * T * X + k3 * D + kmenos1 * (C0 + C1 + C2);
    dXdt = -k1 * T * X + w * C2 + kmenos1 * (C0 + C1 + C2) - k2 * Tp * X + kmenos2 * D + k3 * D;
    dTpdt = w * C2 - k2 * Tp * X + kmenos2 * D;
    dDdt = k2 * Tp * X - (kmenos2 + k3) * D;
    dC0dt = k1 * T * X - (kmenos1 + w) * C0;
    dC1dt = -kmenos1 * C1 - w * C1 + w * C0;
    dC2dt = -kmenos1 * C2 - w * C2 + w * C1;

    % Vector de derivadas
    dydt = [dTdt; dXdt; dTpdt; dDdt; dC0dt; dC1dt; dC2dt];
end


% Función que define el sistema de ecuaciones diferenciales
function dydt = mck(t, y, kon, koff, kp)
    P = y(1);
    T = y(2);
    C0 = y(3);
    C1 = y(4);
    C2 = y(5);
    C3 = y(6);
    C4 = y(7);

    % Ecuaciones diferenciales
    dPdt = -kon * P * T + koff * (C0 + C1 + C2 + C3 + C4);
    dTdt = -kon * P * T + koff * (C0 + C1 + C2 + C3 + C4);
    dC0dt = kon * P * T - (koff + kp) * C0;
    dC1dt = kp * C0 - (koff + kp) * C1;
    dC2dt = kp * C1 - (koff + kp) * C2;
    dC3dt = kp * C2 - (koff + kp) * C3;
    dC4dt = kp * C3 - koff * C4;

    % Vector de derivadas
    dydt = [dPdt; dTdt; dC0dt; dC1dt; dC2dt; dC3dt; dC4dt];
end
