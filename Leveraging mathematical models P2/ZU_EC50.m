%--------------------------------------------------------------------------
% Zero Ultraspecificity model.
% The model is taken from:
%--------------------------------------------------------------------------
% Goldbeter, A., & Koshland Jr, D. E. (1981). An amplified sensitivity arising from covalent modification in biological systems. 
% Proceedings of the National Academy of Sciences, 78(11), 6840-6844.
%--------------------------------------------------------------------------
clear all;

syms kon koff kplus km lambda kp
p = [kon; koff; kplus; km; lambda; kp];

syms T L Tp D C0 P
x = [T; L; Tp; D; C0; P];

f = [
    -kon * T * L + lambda * D + koff * C0;
    -kon * T * L + kp * C0 + koff * C0;
    kon * T* L - (koff + kp) * C0;
    kplus * Tp * P - (km + lambda) * D;
    -kplus * Tp * P + km * D + kp * C0;
    -kplus * Tp * P + km * D + lambda * D
];

%% Emax

final_solution_copia
coeff1_1
coeff1_2
coeff1_3
coeff1_4
coeff1_5
coeff2_1
coeff2_2
coeff2_3
coeff2_4

c11 = coeff1_1var;
c12 = coeff1_2var;
c13 = coeff1_3var;
c14 = coeff1_4var;
c15 = coeff1_5var;
c21 = coeff2_1var;
c22 = coeff2_2var;
c23 = coeff2_3var;
c24 = coeff2_4var;

h = [exprSinDs];



%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0,0,0]; 

u = [];
w = [];
save('ZUEmax','x','p','u','w','h','f','ics','known_ics');

