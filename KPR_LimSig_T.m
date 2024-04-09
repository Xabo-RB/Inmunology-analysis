%--------------------------------------------------------------------------
% KPR McKeithan model.
% The model is taken from:
%--------------------------------------------------------------------------
% Lever, M., Maini, P. K., Van Der Merwe, P. A., & Dushek, O. (2014). 
% Phenotypic models of T cell activation. Nature Reviews Immunology, _14_(9), 619-629.
%--------------------------------------------------------------------------
clear all;


% 4 parameters:
syms kp koff kon phi
p = [kp koff kon phi].';

%% Different steps on KPR

% syms P T C0 C1 C2
% x = [P T C0 C1 C2].';
% N = 1;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * (T+C0+C1+C2); T];
% f = [ 
% 	-kon * P * T + koff * C0 + koff * C1 + koff * C2;
%     -kon * P * T + koff * C0 + koff * C1 + koff * C2;
%     kon * P * T - (koff + kp) * C0;
%     kp * C0 - (koff + phi) * C1;
%     phi * C1 - koff * C2 
% ];

% % N = 2
% syms P T C0 C1 C2 C3
% x = [P T C0 C1 C2 C3].';
% N = 2;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3)); T];
% f = [ 
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3;
%     - kon * P * T + koff*C0 + koff*C1 + koff*C2 + koff*C3;
%     kon * P * T - (koff + kp)*C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + phi)*C2;
%     phi*C2 - (koff)*C3
% ];
% 
% % N = 3
% syms P T C0 C1 C2 C3 C4
% x = [P T C0 C1 C2 C3 C4].';
% N = 3;
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4)); T];
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
% h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4+C5)); T];
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
syms P T C0 C1 C2 C3 C4 C5 C6
x = [P T C0 C1 C2 C3 C4 C5 C6].';
N = 5;
h = [(koff/(phi+koff))*((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4+C5+C6)); T];
f = [
    C0*koff + C1*koff + C2*koff + C3*koff + C4*koff + C5*koff + C6*koff - P*T*kon;
    -kon*P*T + koff*C0 + koff*C1 + koff*C2 + koff*C3 + koff*C4 + koff*C5 + koff*C6; 
    kon * P * T - (koff + kp)*C0; 
    kp*C0 - (koff + kp)*C1; 
    kp*C1 - (koff + kp)*C2; 
    kp*C2 - (koff + kp)*C3; 
    kp*C3 - (koff + kp)*C4; 
    kp*C4 - (koff + phi)*C5; 
    phi*C5 - (koff)*C6
    ];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0]; 

u = [];
w = [];
save('KPRLimSig','x','p','u','w','h','f','ics','known_ics');