%--------------------------------------------------------------------------
% KPR McKeithan model.
% The model is taken from:
%--------------------------------------------------------------------------
% # McKeithan, T. W. (1995). Kinetic proofreading in T-cell receptor signal transduction. 
% #_Proceedings of the national academy of sciences_, _92_(11), 5042-5046.
%--------------------------------------------------------------------------
clear all;


% 4 parameters:
syms kp koff kon
p = [kp koff kon].';

%% Different steps on KPR
% 
syms P T C0 C1
x = [P T C0 C1].';
N = 1;
h = [((kp/(kp+koff))^N) * (T+C0+C1)];
f = [ 
	-kon * P * T + koff * C0 + koff * C1;
    -kon * P * T + koff * C0 + koff * C1;
    kon * P * T - (koff + kp) * C0;
    kp * C0 - koff * C1
];

% % N = 2
% syms P T C0 C1 C2
% x = [P T C0 C1 C2].';
% N = 2;
% h = [((kp/(kp+koff))^N) * ((T+C0+C1+C2))];
% f = [ 
% 	-kon * P * T + koff * C0 + koff * C1 + koff * C2;
%     -kon * P * T + koff * C0 + koff * C1 + koff * C2;
%     kon * P * T - (koff + kp) * C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - koff*C2
% ];
% 
% % N = 3
% syms P T C0 C1 C2 C3
% x = [P T C0 C1 C2 C3].';
% N = 3;
% h = [((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3))];
% f = [ 
% 	-kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3;
%     -kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3;
%     kon * P * T - (koff + kp) * C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff)*C3
% ];
% 
% % N = 4
% syms P T C0 C1 C2 C3 C4
% x = [P T C0 C1 C2 C3 C4].';
% N = 4;
% h = [((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4))];
% f = [ 
% 	-kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3 + koff * C4;
%     -kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3 + koff * C4;
%     kon * P * T - (koff + kp) * C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff + kp)*C3;
%     kp*C3- (koff)*C4
% ];
%
% syms P T C0 C1 C2 C3 C4 C5
% x = [P T C0 C1 C2 C3 C4 C5].';
% N = 5;
% h = [((kp/(kp+koff))^N) * ((T+C0+C1+C2+C3+C4+C5))];
% f = [ 
% 	-kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3 + koff * C4 + koff * C5;
%     -kon * P * T + koff * C0 + koff * C1 + koff * C2 + koff * C3 + koff * C4 + koff * C5;
%     kon * P * T - (koff + kp) * C0;
%     kp*C0 - (koff + kp)*C1;
%     kp*C1 - (koff + kp)*C2;
%     kp*C2 - (koff + kp)*C3;
%     kp*C3 - (koff + kp)*C4;
%     kp*C4 - (koff)*C5
% ];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0]; 

u = [];
w = [];
save('KPRmck','x','p','u','w','h','f','ics','known_ics');
