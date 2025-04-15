%--------------------------------------------------------------------------
% KPR Concentration compensation model.
% The model is taken from:
%--------------------------------------------------------------------------
% Kajita, M. K., Aihara, K., & Kobayashi, T. J. (2017). Balancing specificity, sensitivity, and speed of ligand 
% discrimination by zero-order ultraspecificity. Physical Review E, 96(1), 012405.
%--------------------------------------------------------------------------
clear all;

syms k1 kmenos1 k2 kmenos2 k3 w
p = [k1; kmenos1; k2; kmenos2; k3; w];

syms T X Tp D C0 C1 C2
x = [T; X; Tp; D; C0; C1; C2];

f = [
    -k1 * T * X + k3 * D + kmenos1 * C0 + kmenos1 * C1 + kmenos1 * C2;
    -k1 * T * X + w * C2 + kmenos1 * C0 + kmenos1 * C1 + kmenos1 * C2 - k2 * Tp * X + kmenos2 * D + k3 * D;
    w * C2 - k2 * Tp * X + kmenos2 * D;
    k2 * Tp * X - (kmenos2 + k3) * D;
    k1 * T * X - (kmenos1 + w) * C0;
    -kmenos1 * C1 - w * C1 + w * C0;
    -kmenos1 * C2 - w * C2 + w * C1
];

%% EC50
% EC50
h = [(-(C0 + C1 + C2) * (w * (w / (kmenos1 + w))^2 * ...
      (4 * (C0 + C1 + C2) * k2 * k3 - k2 * k3 * (kmenos1 + 4 * k1 * (T + Tp + C0 + C1 + C2 +D)) + kmenos1 * k1 * (k3 + kmenos2)) + ...
      2 * k2 * k3 * (kmenos1 + w) * (-2 * (C0 + C1 + C2) + kmenos1 + 2 * k1 * (T + Tp + C0 + C1 + C2 +D)))) / ...
      (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 + D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w))...
    + (-(C0 + C1 + C2) * (k2 * k3 * w - k2 * k3 * w * (w / (kmenos1 + w))^2 + k2 * k3 * kmenos1) * ...
      sqrt((kmenos1 * w * (w / (kmenos1 + w))^2 * (k2 * k3 - k1 * (k3 + kmenos2)) - 2 * k2 * k3 * kmenos1 * (kmenos1 + w))^2 / ...
      (k2^2 * k3^2 * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w)^2))) / ...
      (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 +D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w))];

% %% EC50 y TT
% 
% h = [(-(C0 + C1 + C2) * (w * (w / (kmenos1 + w))^2 * ...
%       (4 * (C0 + C1 + C2) * k2 * k3 - k2 * k3 * (kmenos1 + 4 * k1 * (T + Tp + C0 + C1 + C2 +D)) + kmenos1 * k1 * (k3 + kmenos2)) + ...
%       2 * k2 * k3 * (kmenos1 + w) * (-2 * (C0 + C1 + C2) + kmenos1 + 2 * k1 * (T + Tp + C0 + C1 + C2 +D)))) / ...
%       (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 +D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w))...
%     + (-(C0 + C1 + C2) * (k2 * k3 * w - k2 * k3 * w * (w / (kmenos1 + w))^2 + k2 * k3 * kmenos1) * ...
%       sqrt((kmenos1 * w * (w / (kmenos1 + w))^2 * (k2 * k3 - k1 * (k3 + kmenos2)) - 2 * k2 * k3 * kmenos1 * (kmenos1 + w))^2 / ...
%       (k2^2 * k3^2 * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w)^2))) / ...
%       (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 +D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w));
%       (T + Tp + C0 + C1 + C2 +D)];
% 
% %% EC50 y T
% 
% h = [(-(C0 + C1 + C2) * (w * (w / (kmenos1 + w))^2 * ...
%       (4 * (C0 + C1 + C2) * k2 * k3 - k2 * k3 * (kmenos1 + 4 * k1 * (T + Tp + C0 + C1 + C2 +D)) + kmenos1 * k1 * (k3 + kmenos2)) + ...
%       2 * k2 * k3 * (kmenos1 + w) * (-2 * (C0 + C1 + C2) + kmenos1 + 2 * k1 * (T + Tp + C0 + C1 + C2 +D)))) / ...
%       (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 +D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w))...
%     + (-(C0 + C1 + C2) * (k2 * k3 * w - k2 * k3 * w * (w / (kmenos1 + w))^2 + k2 * k3 * kmenos1) * ...
%       sqrt((kmenos1 * w * (w / (kmenos1 + w))^2 * (k2 * k3 - k1 * (k3 + kmenos2)) - 2 * k2 * k3 * kmenos1 * (kmenos1 + w))^2 / ...
%       (k2^2 * k3^2 * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w)^2))) / ...
%       (4 * k2 * k3 * ((C0 + C1 + C2) - k1 * (T + Tp + C0 + C1 + C2 +D)) * (-w * (w / (kmenos1 + w))^2 + kmenos1 + w));
%       T];


%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0,0,0,0]; 

u = [];
w = [];
save('KPCEC50','x','p','u','w','h','f','ics','known_ics');

