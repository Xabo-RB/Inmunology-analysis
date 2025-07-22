%--------------------------------------------------------------------------
% KPR Negative Feedback
% The model is taken from:
%--------------------------------------------------------------------------
% François, P., Voisinne, G., Siggia, E. D., Altan-Bonnet, G., & Vergassola, M. (2013). 
% Phenotypic model for early T-cell activation displaying sensitivity, specificity, and antagonism. 
% Proceedings of the National Academy of Sciences, 110(10), E888-E897.
%--------------------------------------------------------------------------
clear all;


syms kon koff kp b gamma alpha ST beta
p =[kon; koff; kp; b; gamma; alpha; ST; beta];


syms P T C0 C1 C2 S
x = [P T C0 C1 C2 S].';

%% --------------------------------------------- Emax ---------------------------------------------
% q + 
q_plus  = (kp + b + gamma*P + koff + sqrt((kp + b + gamma*P + koff)^2 - 4*kp*(b + gamma*P)))/(2*(b + gamma*P));

% q − 
q_minus = (kp + b + gamma*P + koff - sqrt((kp + b + gamma*P + koff)^2 - 4*kp*(b + gamma*P)))/(2*(b + gamma*P));

h = [(1 - q_minus/q_plus) / (1 - (q_minus/q_plus)^3)* q_minus^2*(T + C0 + C1 + C2)];


%% --------------------------------------------- Ecuaciones ---------------------------------------------
f = [-kon*P*T + koff*(C0 + C1 + C2);                          % dP/dt
    -kon*P*T + koff*(C0 + C1 + C2);                          % dT/dt
     kon*P*T - (koff + kp)*C0 + (b + gamma*S)*C1;                             % dC0/dt
     kp*C0  - (koff + kp + b + gamma*S)*C1 + (b + gamma*S)*C2;                     % dC1/dt
     kp*C1  - (koff + b + gamma*S)*C2;                                   % dC2/dt
     alpha*C1*(ST - S) - beta*S                                                     % dS/dt
];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0]; 

u = [];
w = [];
save('Negative_EC50','x','p','u','w','h','f','ics','known_ics');