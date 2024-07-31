function solution = sensitivity(x0, p, d, tspan)

    neg = @(t,y)negativeII(t, x, p);
    options = odeset('RelTol',1e-3,'AbsTol',1e-3);
    [t,x] = ode45(neg, tspan, x0, options);
    
    lp = length(p); ls = length(sol); lx = length(x0);
    solution = cell(1, lx);
    for i = 1:lx
        solution{i} = zeros(ls, lp + 1);
    end

    for j = 1:lx
        solution{j}(:, 1) = x(:, j);
    end
    for j = 1:lp
        p(j) = p(j) + d * im; % Perturba el parámetro
        neg = @(t,y)negativeII(t, x, p);
        options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        [t,x] = ode45(neg, tspan, x0, options);
        p(j) = complex(real(p(j)), 0);
        x = imag(x) ./ d;

    end



end

function dx = negativeII(t, x, p)
    dx = zeros(7,1);

    dx(1) = p(1) * (p(2) - x(1) - x(2) - x(3)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(5)) + p(4)) * x(1) + (p(7) + p(8)*x(4))*x(2);
    dx(2) = p(4)*x(1) + (p(7) + p(8)*x(4))*x(3) - ((1/p(5)) + p(4) + p(7) + p(8)*x(4))*x(2);
    dx(3) = p(4)*x(2) - ((1/p(5)) + p(7) + p(8)*x(4))*x(3);
    dx(4) = p(10)*(x(2) + x(6))*(p(9) - x(4)) - p(11)*x(4);
    dx(5) = p(1) * (p(12) - x(5) - x(6) - x(7)) * (p(3) - x(1) - x(2) - x(3) - x(5) - x(6) - x(7)) - ((1/p(6)) + p(4)) * x(5) + (p(7) + p(8)*x(4))*x(6);
    dx(6) = p(4)*x(5) + (p(7) + p(8)*x(4))*x(7) - ((1/p(6)) + p(4) + p(7) + p(8)*x(4))*x(6);
    dx(7) = p(4)*x(6) - ((1/p(6)) + p(7) + p(8)*x(4))*x(7);

end