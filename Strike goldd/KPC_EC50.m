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
h = [(-(C0 + C1 + C2) * (w * (w / (koff + w))^N * ...
      (4 * (C0 + C1 + C2) * k2 * k3 - k2 * k3 * (koff + 4 * kon * (T + Tp + C0 + C1 + C2)) + koff * kon * (k3 + km2)) + ...
      2 * k2 * k3 * (koff + w) * (-2 * (C0 + C1 + C2) + koff + 2 * kon * TT))) / ...
      (4 * k2 * k3 * ((C0 + C1 + C2) - kon * TT) * (-w * (w / (koff + w))^N + koff + w))...
    + (-(C0 + C1 + C2) * (k2 * k3 * w - k2 * k3 * w * (w / (koff + w))^N + k2 * k3 * koff) * ...
      sqrt((koff * w * (w / (koff + w))^N * (k2 * k3 - kon * (k3 + km2)) - 2 * k2 * k3 * koff * (koff + w))^2 / ...
      (k2^2 * k3^2 * (-w * (w / (koff + w))^N + koff + w)^2))) / ...
      (4 * k2 * k3 * ((C0 + C1 + C2) - kon * TT) * (-w * (w / (koff + w))^N + koff + w))];
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
known_ics = [0,0,0,0,0,0,0]; 

u = [];
w = [];
save('KPCEC50','x','p','u','w','h','f','ics','known_ics');

