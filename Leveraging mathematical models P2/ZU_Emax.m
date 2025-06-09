%--------------------------------------------------------------------------
% Zero Ultraspecificity model.
% The model is taken from:
%--------------------------------------------------------------------------
% Goldbeter, A., & Koshland Jr, D. E. (1981). An amplified sensitivity arising from covalent modification in biological systems. 
% Proceedings of the National Academy of Sciences, 78(11), 6840-6844.
%--------------------------------------------------------------------------
clear all;

syms k1 kmenos1 k2 kmenos2 k3 w
p = [k1; kmenos1; k2; kmenos2; k3; w];

syms T L Tp D C0 P
x = [T; L; Tp; D; C0; P];

f = [
    -k1 * T * L + k3 * D + kmenos1 * C0;
    -k1 * T * L + w * C0 + kmenos1 * C0;
    k1 * T* L - (kmenos1 + w) * C0;
    k2 * Tp * P - (kmenos2 + k3) * D;
    -k2 * Tp * P + kmenos2 * D + w * C0;
    -k2 * Tp * P + kmenos2 * D + k3 * D
];

%% Emax

h = [( (w*k2*k3*(T + Tp + C0 + D) - k2*k3*(k3 + w)*(D + P) - w*k3*(kmenos2 + k3)) + ...
    sqrt( (w*k2*k3*(T + Tp + C0 + D) - k2*k3*(k3 + w)*(D + P) - w*k3*(kmenos2 + k3))^2 ...
             + 4*w^2*k3^2*k2*(kmenos2 + k3)*(T + Tp + C0 + D) ) )/ ( 2*w*k2*k3 ); T + Tp + C0 + D];



%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0,0,0]; 

u = [];
w = [];
save('ZUEmax','x','p','u','w','h','f','ics','known_ics');

