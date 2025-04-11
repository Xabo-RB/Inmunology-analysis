%--------------------------------------------------------------------------
% KPR Sustained signaling model.
% The model is taken from:
%--------------------------------------------------------------------------
% Coombs, D., Kalergis, A. M., Nathenson, S. G., Wofsy, C., & Goldstein, B. (2002).
% Activated TCRs remain marked for internalization after dissociation from pMHC.
% Nature Immunology, 3(10), 926-931.

% González, P. A., Carreño, L. J., Coombs, D., Mora, J. E., Palmieri, E., Goldstein, B., ... & Kalergis, A. M. (2005).
% T cell receptor binding kinetics required for T cell activation depend on the density of cognate ligand on the antigen-presenting cell.
% Proceedings of the National Academy of Sciences, 102(13), 4824-4829.
%--------------------------------------------------------------------------
clear all;


% 4 parameters:
syms kp koff kon lambda
p = [kp koff kon lambda].';

%% Different steps on KPR

syms P T C0 C1 C2 Tast
x = [P T C0 C1 C2 Tast].';
N = 2;
h = [ (((kp/(kp + koff))^N) * (C0+C1+C2) * (koff + lambda - kon*(C0+C1+C2) + kon*0.5*(T+C0+C1+C2+Tast)) - (lambda * (T+C0+C1+C2+Tast) *0.5)) / ( ((kp/(kp + koff))^N) * kon * ((T+C0+C1+C2+Tast) * 0.5 - (C0+C1+C2)) ) ];
f = [ 
	- kon * P * T + koff * C0 + koff * C1 + koff * C2 - kon * P * Tast;
    - kon * P * T + koff * C0 + koff * C1 + lambda * Tast;
    kon * P * T - (koff + kp) * C0;
    kp * C0 - (koff + kp) * C1;
    kp * C1 - koff * C2 + kon * P * Tast;
    koff * C2 - kon * P * Tast - lambda * Tast
];

% 
% % N = 3
% syms P T C0 C1 C2 C3 C4
% x = [P T C0 C1 C2 C3 C4].';
% N = 3;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4))];
% f = [ 
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4;
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4;
%     kon * P * T - (koff + kp)*C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff + phi)*C3;
%     phi*C3 - (koff)*C4
% ];
% 
% % N = 4
% syms P T C0 C1 C2 C3 C4 C5
% x = [P T C0 C1 C2 C3 C4 C5].';
% N = 4;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4+C5))];
% f = [ 
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4 +(koff)*C5;
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4 +(koff)*C5;
%     kon * P * T - (koff + kp)*C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff + kp)*C3;
%     kp*C3 - (koff + phi)*C4;
%     phi*C4 - (koff)*C5
% ];
%
% syms P T C0 C1 C2 C3 C4 C5 C6
% x = [P T C0 C1 C2 C3 C4 C5 C6].';
% N = 5;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4+C5+C6))];
% f = [ 
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4 +(koff)*C5 +(koff)*C6;
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3 +(koff)*C4 +(koff)*C5 +(koff)*C6;
%     kon * P * T - (koff + kp)*C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff + kp)*C3;
%     kp*C3 - (koff + kp)*C4;
%     kp*C4 - (koff + phi)*C5;
%     phi*C5 - (koff)*C6
% ];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0,0,0]; 

u = [];
w = [];
save('KPRSustSig','x','p','u','w','h','f','ics','known_ics');