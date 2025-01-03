syms S T A
x = [S; T; A];

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms lambda phi ki h keff L kon koff L0 sigma
p =[lambda; phi; ki; h; keff; kon; koff; L0; sigma];

% 1 output
%h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% Linfocitos T(t)
% h = T;

%TT
h = T + A + S;

% Response R(t)
% h = A;

% Linfocitos y Respuesta
%h = [T; A];

%TT y Respuesta
%h = [T + A + S; A];

% Emax para h = 1
%h = (s * phi + s * lambda * phi + s*s) / (ki * (lambda * phi + s));


% dynamic equations
% f = [-lambda * phi * (S-T) + s*(1-S);
%     phi * (S - T) + s*(1-T) - k * (T^hh)*(L^hh);
%     k *(T^hh)*(L^hh) - ki*A];

% dynamic equations
f = [-phi*((1 + lambda)/lambda)*(lambda*S - T) + (sigma/(1 + lambda))*(100 - S - T);
    phi*((1 + lambda)/lambda)*(lambda*S - T) + sigma*(lambda/(1 + lambda))*(100 - S - T) - keff*T^h*(L0*koff/(koff + kon*T))^h;
    keff*T^h*(L0*koff/(koff + kon*T))^h - ki*A];

% syms lambda phi ki keff kon koff L0 sigma
% p =[lambda; phi; ki; keff; kon; koff; L0; sigma];
% 
% % dynamic equations
% f = [-phi*((1 + lambda)/lambda)*(lambda*S - T) + (sigma/(1 + lambda))*(100 - S - T);
%     phi*((1 + lambda)/lambda)*(lambda*S - T) + sigma*(lambda/(1 + lambda))*(100 - S - T) - keff*T^2*(L0*koff/(koff + kon*T))^2;
%     keff*T^2*(L0*koff/(koff + kon*T))^2 - ki*A];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('SerialTriggering','x','p','h','f','u','ics','known_ics');


