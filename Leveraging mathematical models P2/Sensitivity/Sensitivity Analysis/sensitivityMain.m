function solution = sensitivityMain(x0, p, d, tspan)

    neg = @(t,y)negativeII(t, y, p);
    options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
    [t,x] = ode45(neg, tspan, x0, options);
    
    lp = length(p); ls = size(x, 1); lx = length(x0);
    solution = cell(1, lx);
    for i = 1:lx
        solution{i} = zeros(ls, lp + 1);
    end

    for j = 1:lx
        solution{j}(:, 1) = x(:, j);
    end

    for j = 1:lp

        p(j) = p(j) + d * 1i;
        
        options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'Refine', 1);
        neg = @(t,y)negativeII(t, y, p);
        [t,x] = ode45(neg, tspan, x0, options);
        
        p(j) = complex(real(p(j)), 0);
        xSens = imag(x) ./ d;
        
        for k = 1:lx
            solution{k}(:, j + 1) = xSens(:, k);
        end


    end

end