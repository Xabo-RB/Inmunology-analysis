%--------------------------------------------------------------------------
% Serial triggering
% The model is taken from:
%--------------------------------------------------------------------------
% Salvatore Valitutti, Sabina M¨uller, Marina Cella, Elisabetta Padovan, and Antonio Lanza-
%vecchia. Serial triggering of many t-cell receptors by a few peptide–mhc complexes. Nature,
%375(6527):148–151, 1995
%--------------------------------------------------------------------------
clear all;


syms lambda phi ki keff L sigma
p =[lambda; phi; ki; keff; L; sigma];


syms S T Tp
x = [S T Tp].';

%% --------------------------------------------- Emax ---------------------------------------------
h = [sigma * ( sigma*lambda + phi*(1 + lambda) )/( ki * (1 + lambda) * (sigma + phi) ); S + T + Tp];

%% --------------------------------------------- Normal ---------------------------------------------
h = [T];


f = [ 
	-phi*(S - T/lambda) + sigma*(1/(1 + lambda) - S);
    phi*(S - T/lambda) + sigma*(lambda/(1 + lambda) - T) - keff*T^2*L^2;
    keff*T^2*L^2 - ki*Tp;
];

%% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0]; 

u = [];
w = [];
save('SerialTrig_Emax','x','p','u','w','h','f','ics','known_ics');