k_plus = 1;
lambda = 1;
k_minus = 1;
k_p = 1;
X_T = 1;
T_T = 1;
kon = 1;
koff = 1;
N = 5;
psi = k_p / (koff + k_p);
Phi = (k_p / lambda) * ((psi^N - psi^(N + 1)) / (1 - psi^(N + 1)));


beta1 = k_plus * lambda * Phi + ((k_minus + lambda) * Phi + k_plus * (X_T - T_T) * (1 + Phi)) * kon * psi^(N + 1);
beta2 = k_plus * lambda * Phi + ((k_minus + lambda) * Phi + k_plus * (X_T + T_T) * (1 + Phi)) * kon * psi^(N + 1);


numerador = -(k_minus + lambda) * (beta1 - sqrt(beta2^2 - 4 * k_plus^2 * kon^2 * T_T * X_T * (1 + Phi)^2 * psi^(2 * (N + 1))));
denominador = 2 * k_plus * (1 + Phi) * (k_plus * lambda + kon * (k_minus + lambda)) * psi^(N + 1);
resultado = numerador / denominador;

b1 = k_plus * (Phi * lambda + koff) + (Phi + 1) * k_plus * kon * (T_T + X_T) + Phi * kon * (k_minus + lambda);
b2 = k_plus * (Phi * lambda + koff) + Phi * (k_minus + lambda) * kon;
Tp_hat = Phi * (k_minus + lambda) * ((-b1 + sqrt(b1^2 - 4 * ((Phi + 1) * k_plus * kon)^2 * T_T * X_T)) / ...
    (2 * b2 * (1 + Phi) * k_plus));
