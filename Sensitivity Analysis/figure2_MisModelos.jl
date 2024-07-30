using Statistics, Random, LinearAlgebra, Plots, QuadGK, Polyester, Enzyme
using DifferentialEquations
using DataFrames, LaTeXStrings, ForwardDiff, DiffEqSensitivity, Sundials, NumericIO
#, DelimitedFiles, #, BenchmarkTools, Revise, PyPlot
# El paquete BioSimulator, AnalyticSensitivity hay que instalarlo a mano desde el repositorio de github
# Pkg.add(url="https://github.com/alanderos91/BioSimulator.jl.git")
# Pkg.add(url="https://github.com/rachelmester/SensitivityAnalysis")
using BioSimulator
using AnalyticSensitivity
using CairoMakie

ASPkg = AnalyticSensitivity
versioninfo()

include("Modelos.jl")

#   OCCUPANCY       (1)
#   McKeithan       (2)
#   McKeithan10     (3)
#   Limited         (4)
#   Sustained       (5)
#   Negative        (6)
#   Galvez          (7)
#   LimIFF          (8)
#   Induced         (9)
#   ST              (10)


case = 5

if case == 1

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEOccupancy, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEOccupancy, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    kon = 5e-5; koff = 0.01; afinidad = kon/koff
    p = complex([5e-5, 0.01]); # kon koff  (afinidad = 0.005)
    x0 = complex([0, 100, 2e4]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    solution = sensitivity(x0, p, d, tspan); # find solution and partials

    #==
    (d, tspan) = (1.0e-16, (0.0,100)); # step size and time interval in days
    p1 = complex([0.0005, 0.1]); # kon koff (afinidad = 0.005)
    solution1 = sensitivity(x0, p1, d, tspan); # find solution and partials

    (d, tspan) = (1.0e-16, (0.0,10000)); # step size and time interval in days
    p2 = complex([5e-6, 0.001]); # kon koff (afinidad = 0.005)
    solution2 = sensitivity(x0, p2, d, tspan); # find solution and partials   

    ss1 = log10.(abs.(solution[1][:, 3]))
    ss2 = log10.(abs.(solution1[1][:, 3]))
    ss3 = log10.(abs.(solution2[1][:, 3]))
    pTlog = Plots.plot(ss1, label = "koff = 0.01", xlabel= "t", ylabel = "R(t)")
    Plots.plot!(ss2, label = "koff = 0.1")
    Plots.plot!(ss3, label = "koff = 0.001")
    display(pTlog)

    pT = plot(solution[1][:, 3], label = "koff = 0.01", xlabel= "t", ylabel = "R(t)")
    Plots.plot!(solution1[1][:, 3], label = "koff = 0.1")
    Plots.plot!(solution2[1][:, 3], label = "koff = 0.001")
    #Plots.plot(solution[1][:, 3], label = "x1", xlabel= "t", ylabel = "R(t)") #xlims = (tspan[1],tspan[2]))
    display(pT)

    p1 = plot(solution[1][:, 3], label = "koff = 0.01", xlabel= "t", ylabel = "R(t)", title = "Solución 1")
    display(p1)
    p2 = plot(solution1[1][:, 3], label = "koff = 0.1", xlabel= "t", ylabel = "R(t)", title = "Solución 2")
    display(p2)
    p3 = plot(solution2[1][:, 3], label = "koff = 0.001", xlabel= "t", ylabel = "R(t)", title = "Solución 3")
    display(p3)
    ==#
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR

    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)

        #konCorregido = afinidad*koffVect[i]
        #p = complex([konCorregido, koffVect[i]]);
        p = complex([5e-5, koffVect[i]]);
        solution = sensitivity(x0, p, d, tspan);

        newSol = (solution[1][:, 3].*koffVect[i])./solution[1][:, 1]
        results_matrix[i, :] = newSol

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =50, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        title = "Occupancy model", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig


elseif case == 2 # ==============================================================================================================

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEKPRmcK, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEKPRmcK, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days


    p = complex([5e-5, 0.01, 1]); # kon koff kp
    solution = sensitivity(x0, p, d, tspan); 
    #== PLOTEAR LA TRAYECTORIA ÚNICA DEL OUTPUT
    p1 = Plots.plot(solution[4][:, 3], label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    ==# 

    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[4][:, 3]))
    for i in eachindex(koffVect)

        p = complex([5e-5, koffVect[i], 1]);
        solution = sensitivity(x0, p, d, tspan);

        newSol = (solution[4][:, 3].*koffVect[i])./solution[4][:, 1]
        results_matrix[i, :] = newSol

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =50, step = 1))

    #== HEATMAP con plots sin interpolar
    p = Plots.plot(Plots.heatmap(time1, koffVect, results_matrix),
    xlabel="Time (s)", 
    ylabel= "Dissociate rate",
    title= L"Response sensitivity to $k_{off}$",
    interpolate=true)
    display(p)
    ==#

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        title = "McKeithan", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig


elseif case == 3 # ==============================================================================================================


    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEKPRmcK10, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEKPRmcK10, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    p = complex([5e-5, 0.01, 1]); # kon koff kp
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]); # initial values
    #x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    solution = sensitivity(x0, p, d, tspan); # find solution and partials
    Plots.plot(solution[13][:, 1], label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))

elseif case == 4

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEKPRLimSig, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEKPRLimSig, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,100)); # step size and time interval in days
    p = complex([5e-5, 0.01, 1, 0.09]); # kon koff kp phi
    solution = sensitivity(x0, p, d, tspan); 
    #==
    # LA SALIDA ES C_N, es decir, x8
    p1 = Plots.plot(solution[8][:, 1], label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    ==#
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[8][:, 3]))
    for i in eachindex(koffVect)

        p = complex([5e-5, koffVect[i], 1, 0.09]);
        solution = sensitivity(x0, p, d, tspan);

        newSol = (solution[8][:, 3].*koffVect[i])./solution[8][:, 1]
        results_matrix[i, :] = newSol

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        title = "Limited signaling", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    

elseif case == 5

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEKPRSustSig, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEKPRSustSig, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    #kon = 5e-5; koff = 0.01; afinidad = kon/koff
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,10000)); # step size and time interval in days
    p = complex([5e-5, 0.01, 1, 0.001]); # kon koff kp lambda
    solution = sensitivity(x0, p, d, tspan); 
    
    #==
    NewSolR = solution[8][:, 1] .+ solution[9][:, 1]
    # LA SALIDA ES C_N + Tast, es decir, x8 + x9
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    ==#
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[8][:, 3]))
    for i in eachindex(koffVect)

        #konCorreg = afinidad*koffVect[i]
        #p = complex([konCorreg, koffVect[i], 1, 0.001]);
        #p = complex([5e-5, koffVect[i], 1, 0.001]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[8][:, 3] .+ solution[9][:, 3]
        SolResponseEst = solution[8][:, 1] .+ solution[9][:, 1]
        newSol = (SolResponse.*koffVect[i])./SolResponseEst
        results_matrix[i, :] = newSol

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =10000, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Sustained signaling", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 6

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEKPRNegFeed, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEKPRNegFeed, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,200)); # step size and time interval in days
    p = complex([5e-5, 0.01, 1, 4.4e-4, 0.04, 1, 2e-4, 6e5]); # kon koff kp gama b beta alpha ST
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[8][:, 3]
    # LA SALIDA ES C_N + Tast, es decir, x8 + x9
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[8][:, 3]))
    for i in eachindex(koffVect)

        p = complex([5e-5, koffVect[i], 1, 4.4e-4, 0.04, 1, 2e-4, 6e5]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[8][:, 3]
        newSol = (SolResponse.*koffVect[i])./solution[8][:, 1]
        results_matrix[i, :] = newSol


        #results_matrix[i, :] = log10.(abs.(newSol))

        #results_matrix[i, :] = solution[4][:, 3]

    end
    #==
    results_matrix = results_matrix[:, 2:end]
    minimo = abs(minimum(results_matrix))
    results_matrix = minimo .+ results_matrix
    results_matrix = log10.(results_matrix)
    ==#

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =200, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Negative regulator", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 7

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEStabChain2, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEStabChain2, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,100)); # step size and time interval in days
    p = complex([5e-5, 0.01, 1, 1.5, 1.03]); # kon koff kp r rp
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[5][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([5e-5, koffVect[i], 1, 1.5, 1.03]); # kon koff kp r rp
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[5][:, 3]
        newSol = (SolResponse.*koffVect[i])./solution[5][:, 1]
        results_matrix[i, :] = newSol
        #results_matrix[i, :] = log10.(abs.(newSol))

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Stabilizing chain", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 8

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODELimIFF, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODELimIFF, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,100)); # step size and time interval in days
    # phi = p[4],   gamma = p[5],   lambda = p[6],  delta = p[7],   YT = p[8],  PT = p[9],  mu = p[10]
    p = complex([5e-5, 0.01, 1, 0.09, 1, 0.5, 50, 100, 100, 2.5, 500]); # Valores sacados de Galvez
    #p = complex([5e-5, 0.01, exp(-2.6), exp(-0.4), 500, exp(4.4), exp(3.4), 500, 500, exp(4.2), 500]); # Valores sacados del propio paper
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[7][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([5e-5, koffVect[i], 1, 0.09, 1, 0.5, 50, 100, 100, 2.5, 500]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[7][:, 3]
        newSol = (SolResponse.*koffVect[i])./solution[7][:, 1]
        results_matrix[i, :] = newSol
        #results_matrix[i, :] = log10.(abs.(newSol))

        #results_matrix[i, :] = solution[4][:, 3]

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Limited signaling IFF", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 9

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEIndReb, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEIndReb, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    
    x0 = complex([100, 2e4, 0, 0, 0, 0, 0, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    # rho1 = p[4],  lambdaR = p[5], rho2 = p[6]
    p = complex([5e-5, 0.01, 1, 0.09, 1000, 10000, 1000]); # Valores sacados de Galvez
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[8][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([5e-5, koffVect[i], 1, 0.09, 1000, 10000, 1000]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[8][:, 3]
        newSol = (SolResponse.*koffVect[i])./solution[8][:, 1]

        #SolResponse = solution[7][:, 3] + solution[5][:, 3]
        #newSol = (SolResponse.*koffVect[i])./(solution[7][:, 1] + solution[5][:, 1])

        
        results_matrix[i, :] = newSol

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Induced rebinding", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 10

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEST, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEST, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    # S(t) = x1, T(t) = x2, A(t) = x3, Y(t) = x4,
    x0 = complex([0.01, 0, 0, 100]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    # lambda = p[1], phi = p[2], s = p[3], keff = p[4], ki = p[5], L = p[6], h = p[7]
    p = complex([0.61, 0.0055, 0.0011, 1, 0.094, 100, 5]); # Valores sacados de Galvez
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[4][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([0.61, 0.0055, 0.0011, koffVect[i], 0.094, 100, 5]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[4][:, 5]
        newSol = (SolResponse.*koffVect[i])./solution[4][:, 5]

        #SolResponse = solution[7][:, 3] + solution[5][:, 3]
        #newSol = (SolResponse.*koffVect[i])./(solution[7][:, 1] + solution[5][:, 1])

        
        results_matrix[i, :] = newSol

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Induced rebinding", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 11

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODEBachmann, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODEBachmann, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    # S(t) = x1, T(t) = x2, A(t) = x3, Y(t) = x4,
    x0 = complex([100]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    # keff = p[1], . = p[2], . = p[3]
    p = complex([0.6, 100, 0.1]); # Valores sacados de Galvez
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[4][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([0.6, 100, koffVect[i]]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[1][:, 4]
        newSol = (SolResponse.*koffVect[i])./solution[1][:, 1]

        #SolResponse = solution[7][:, 3] + solution[5][:, 3]
        #newSol = (SolResponse.*koffVect[i])./(solution[7][:, 1] + solution[5][:, 1])

        
        results_matrix[i, :] = newSol

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Induced rebinding", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
elseif case == 12

    function sensitivity(x0, p, d, tspan)
        problem = ODEProblem{true}(ODENegativeII, x0, tspan, p)
        sol = solve(problem, saveat = 1.0) # solve ODE
        (lp, ls, lx) = (length(p), length(sol), length(x0))  
        solution = Dict{Int, Any}(i => zeros(ls, lp + 1) for i in 1:lx)
        for j = 1:lx # record solution for each species
            @views solution[j][:, 1] = sol[j, :]
        end
        for j = 1:lp
            p[j] = p[j] + d * im # perturb parameter
            problem = ODEProblem{true}(ODENegativeII, x0, tspan, p)
            sol = solve(problem, saveat = 1.0) # resolve ODE
            p[j] = complex(real(p[j]), 0.0) # reset parameter
            @views sol .= imag(sol) / d # compute partial
            for k = 1:lx # record partial for each species
            @views solution[k][:,j + 1] = sol[k, :]
            end
        end
        return solution
    end

    # S(t) = x1, T(t) = x2, A(t) = x3, Y(t) = x4,
    x0 = complex([100]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    # keff = p[1], . = p[2], . = p[3]
    p = complex([0.6, 100, 0.1]); # Valores sacados de Galvez
    solution = sensitivity(x0, p, d, tspan); 
    
    
    NewSolR = solution[4][:, 1]
    p1 = Plots.plot(NewSolR, label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))
    display(p1)
    
    
    # ------------- RECOGER LOS RESULTADOS DE SENSIBILIDAD PARA CADA KOFF DEL VECTOR
    koffVect = collect(range(0.001, stop =1, step = 0.001))
    results_matrix = zeros(length(koffVect), length(solution[1][:, 3]))
    for i in eachindex(koffVect)
        p = complex([0.6, 100, koffVect[i]]);
        solution = sensitivity(x0, p, d, tspan);

        SolResponse = solution[1][:, 4]
        newSol = (SolResponse.*koffVect[i])./solution[1][:, 1]

        #SolResponse = solution[7][:, 3] + solution[5][:, 3]
        #newSol = (SolResponse.*koffVect[i])./(solution[7][:, 1] + solution[5][:, 1])

        
        results_matrix[i, :] = newSol

    end

    # VECTOR TIEMPO PARA PLOTEAR
    time1 = collect(range(0, stop =100, step = 1))

    # Suponiendo que time1, koffVect y results_matrix ya están definidos
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], 
        #title = L"Response sensitivity to $k_{off}$", 
        title = "Induced rebinding", 
        xlabel = "Time (s)", 
        ylabel = "Dissociate rate"
        )
    hm = CairoMakie.heatmap!(ax, time1, koffVect, results_matrix', interpolate = true, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Sensitivity") 
    fig
    
    
end


