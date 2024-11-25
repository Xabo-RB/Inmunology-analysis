%--------------------------------------------------------------------------
% KPR Concentration compensation model.
% The model is taken from:
%--------------------------------------------------------------------------
% Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of ligand 
% discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.
%--------------------------------------------------------------------------
clear all;

syms k1 kmenos1 k2 kmenos2 k3 w
params = [k1; kmenos1; k2; kmenos2; k3; w];

syms T X Tp D C0 C1 C2
x = [T; X; Tp; D; C0; C1; C2];
% EC50
h = [];
f = [
    -k1 * T * X + k3 * D + kmenos1 * C0 + kmenos1 * C1 + kmenos1 * C2;
    -k1 * T * X + w * C2 + kmenos1 * C0 + kmenos1 * C1 + kmenos1 * C2 - k2 * Tp * X + kmenos2 * D + k3 * D;
    w * C2 - k2 * Tp * X + kmenos2 * D;
    k2 * Tp * X - (kmenos2 + k3) * D;
    k1 * T * X - (kmenos1 + w) * C0;
    -kmenos1 * C1 - w * C1 + w * C0;
    -kmenos1 * C2 - w * C2 + w * C1
];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0]; 

u = [];
w = [];
save('KPRLimSig','x','p','u','w','h','f','ics','known_ics');

