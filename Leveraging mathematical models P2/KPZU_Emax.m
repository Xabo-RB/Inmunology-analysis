%--------------------------------------------------------------------------
% Kinetic Proofreading wit Zero-order Ultraspecificity model.
% The model is taken from:
%--------------------------------------------------------------------------
% Goldbeter, A., & Koshland Jr, D. E. (1981). An amplified sensitivity arising from covalent modification in biological systems. 
% Proceedings of the National Academy of Sciences, 78(11), 6840-6844.
%--------------------------------------------------------------------------
clear all;

syms k1 kmenos1 k2 kmenos2 k3 w
p = [k1; kmenos1; k2; kmenos2; k3; w];

syms T P Tp Q D C0 C1 C2 
x = [T; P; Tp; Q; D; C0; C1; C2];

f = [
    -k1 * T * P + k3 * D + kmenos1 * C0 + kmenos1 * C1 + kmenos1 * C2;
    -k1 * T * P + w*C2 + kmenos1 * C0 + kmenos1*C1 + kmenos1*C2;
    w*C2 - k2*Tp*Q + kmenos2*D;
    -k2 * Tp * Q + kmenos2 * D + k3 * D;
    k2 * Tp * Q - (kmenos2 + k3) * D;
    k1 * T * P - (kmenos1 + w) * C0;
    - kmenos1*C1 - w*C1 + w*C0;
    - kmenos1*C2 - w*C2 + w*C1
];

%% Emax

h = [
    ( ((w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))) *...
    k2 * k3 * (T + Tp + D + C0 + C1 + C2) - k2 * k3 * (k3 + (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))))...
    * (Q + D) - (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))) * k3 * (kmenos2 + k3)) + ...
    sqrt( ((w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))) ...
    * k2 * k3 * (T + Tp + D + C0 + C1 + C2) - k2 * k3 * (k3 + (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3)))) ...
    * (Q + D) - (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))) * k3 ...
    * (kmenos2 + k3))^2 + 4 * (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3)))^2 ...
    * k3^2 * k2 * (kmenos2 + k3) * (T + Tp + D + C0 + C1 + C2) ) ) / ( 2 * (w * (w / (kmenos1 + w))^2 * ((1 - (w / (kmenos1 + w))) / (1 - (w / (kmenos1 + w))^3))) ...
    * k2 * k3 ); T + Tp + D + C0 + C1 + C2
    ];


%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0,0,0,0,0]; 

u = [];
w = [];
save('KPZUEmax','x','p','u','w','h','f','ics','known_ics');
