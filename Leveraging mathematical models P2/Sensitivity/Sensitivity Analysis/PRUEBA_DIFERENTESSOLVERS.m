clear
clc

% === TUS VALORES, sin complex ni d ===
x0 = [100, 2e4, 0, 0, 0, 0, 0, 0, 0, 0]';   % 10x1
tspan = 0.0:0.05:200;

p = [ ...
    5e-5, ...         % p1 (kon)
    0.01, ...         % p2 (koff)
    1, ...            % p3
    4.4e-4, ...       % p4
    0.04, ...         % p5
    1, ...            % p6
    2e-4, ...         % p7
    6e5, ...          % p8
    exp(-2), ...      % p9
    exp(-1), ...      % p10
    exp(0) ...        % p11 = 1
];

% === Tolerancias (recomendado AbsTol vector por escalas distintas) ===
AbsTol = [ ...
    1e-7;  ... L ~ 1e2
    1e-4;  ... P ~ 2e4
    1e-10; ... C0
    1e-10; ... C1
    1e-10; ... C2
    1e-10; ... C3
    1e-10; ... C4
    1e-10; ... C5
    1e-10; ... S
    1e-10  ... Y (estado 10)
];

% === Benchmark (sin Jacobiano) ===
results = benchmarkKPRSolvers_tspan(p, x0, tspan, ...
    'RelTol', 1e-6, 'AbsTol', AbsTol, ...
    'Runs', 3, 'NonNegative', 1:9);

function results = benchmarkKPRSolvers_tspan(p, x0, tspan, varargin)
% BENCHMARKKPRSOLVERS_TSPAN  Compara varios solvers ODE con un tspan dado.
% Uso:
%   results = benchmarkKPRSolvers_tspan(p, x0, tspan)
%   results = benchmarkKPRSolvers_tspan(p, x0, tspan, 'RelTol',1e-6,'AbsTol',1e-9,'Runs',3,'NonNegative',1:9)
%
% Entradas:
%   p     : parámetros (1x11)
%   x0    : condiciones iniciales (10x1)
%   tspan : vector de tiempos (por ej. 0:0.05:200)
%
% Opcionales:
%   'RelTol'      (1e-6 por defecto)
%   'AbsTol'      (1e-9 por defecto; puede ser escalar o vector 10x1)
%   'Runs'        nº repeticiones para mediana de tiempo (3 por defecto)
%   'MaxStep'     paso máximo (default [])
%   'InitialStep' paso inicial (default [])
%   'NonNegative' índices de estados no negativos (default [])
%
% Salida:
%   results: tabla {solver, ok, time_s, fevals, steps, tf_points, error}

    ip = inputParser;
    ip.addParameter('RelTol', 1e-6);
    ip.addParameter('AbsTol', 1e-9);
    ip.addParameter('Runs', 3);
    ip.addParameter('MaxStep', []);
    ip.addParameter('InitialStep', []);
    ip.addParameter('NonNegative', []);
    ip.parse(varargin{:});
    optsIn = ip.Results;

    % Opciones comunes (sin Jacobiano)
    baseOpts = odeset( ...
        'RelTol', optsIn.RelTol, ...
        'AbsTol', optsIn.AbsTol, ...
        'Refine', 1, ...
        'MaxStep', optsIn.MaxStep, ...
        'InitialStep', optsIn.InitialStep ...
    );
    if ~isempty(optsIn.NonNegative)
        baseOpts = odeset(baseOpts, 'NonNegative', optsIn.NonNegative);
    end

    % Lista de solvers
    solvers = {@ode45, @ode23, @ode113, @ode15s, @ode23s, @ode23t, @ode23tb};
    solverNames = cellfun(@func2str, solvers, 'UniformOutput', false);

    resultsCell = cell(numel(solvers), 1);

    for k = 1:numel(solvers)
        solver = solvers{k};
        name   = solverNames{k};

        fevals = 0;
        steps  = 0;

        % ODE envuelta para contar evaluaciones
        odefun_wrapped = @(t,x) count_odefun(t,x);

        % OutputFcn para contar pasos aceptados
        outfun = @(t,y,flag) step_counter(t,y,flag);

        opts = odeset(baseOpts, 'OutputFcn', outfun);

        times = zeros(optsIn.Runs,1);
        last_t = [];
        ok = true; errmsg = '';

        for r = 1:optsIn.Runs
            fevals = 0; steps = 0; last_t = [];
            try
                tic;
                [t,~] = solver(odefun_wrapped, tspan, x0, opts); %#ok<ASGLU>
                times(r) = toc;
                last_t = t;
            catch ME
                ok = false;
                errmsg = ME.message;
                times(r) = NaN;
                break;
            end
        end

        time_s = median(times,'omitnan');
        tf_points = numel(last_t);

        resultsCell{k} = {name, ok, time_s, fevals, steps, tf_points, errmsg};
    end

    results = cell2table(vertcat(resultsCell{:}), ...
        'VariableNames', {'solver','ok','time_s','fevals','steps','tf_points','error'});

    okMask = results.ok & ~isnan(results.time_s);
    %[~,ord] = sort(results.time_s(okMask));
    results = [sortrows(results(okMask,:), 'time_s'); results(~okMask,:)];
    disp('== Benchmark (menor tiempo es mejor) ==');
    disp(results);

    % ======= funciones anidadas =======
    function dx = count_odefun(t,x)
        fevals = fevals + 1;
        dx = ODEKPRNegFeed2(t, x, p);
    end

    function status = step_counter(t,~,flag)
        if isempty(flag)
            steps = steps + 1;
        end
        status = 0; % continuar
    end
end

function dx = ODEKPRNegFeed2(t, x, p)
    % Inicializar el vector dx con ceros
    dx = zeros(10, 1);
    
    % Definir las ecuaciones diferenciales
    dx(1) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8)); % L
    dx(2) = -p(1) * x(1) * x(2) + p(2) * (x(3) + x(4) + x(5) + x(6) + x(7) + x(8)); % P
    dx(3) = p(1) * x(1) * x(2) - (p(2) + p(3)) * x(3) + (p(5) + p(4) * x(9)) * x(4); % C0
    dx(4) = p(3) * x(3) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(4) + (p(5) + p(4) * x(9)) * x(5); % C1
    dx(5) = p(3) * x(4) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(5) + (p(5) + p(4) * x(9)) * x(6); % C2
    dx(6) = p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7); % C3
    dx(7) = p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8); % C4
    dx(8) = p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8); % C5
    dx(9) = p(7) * x(4) * (p(8) - x(9)) - p(6) * x(9); % S
    dx(10) = p(9)*(p(3) * x(5) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(6) + (p(5) + p(4) * x(9)) * x(7)) + ...
        p(10)*(p(3) * x(6) - (p(2) + p(3) + p(5) + p(4) * x(9)) * x(7) + (p(5) + p(4) * x(9)) * x(8)) + ...
        p(11)*(p(3) * x(7) - (p(2) + p(5) + p(4) * x(9)) * x(8));
end

    