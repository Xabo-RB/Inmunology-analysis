using Statistics, Random, LinearAlgebra, Plots, QuadGK, Polyester, Enzyme
using DifferentialEquations
using DataFrames, LaTeXStrings, ForwardDiff, DiffEqSensitivity, Sundials, NumericIO
#, DelimitedFiles, #, BenchmarkTools, Revise, PyPlot
# El paquete BioSimulator, AnalyticSensitivity hay que instalarlo a mano desde el repositorio de github
# Pkg.add(url="https://github.com/alanderos91/BioSimulator.jl.git")
using BioSimulator
using AnalyticSensitivity

#ASPkg = AnalyticSensitivity
#versioninfo()

include("Modelos.jl")

#   OCCUPANCY       (1)
#   McKeithan       (2)
#   McKeithan10     (3)
case = 1

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

    p = complex([5e-5, 0.01]); # kon koff  (afinidad = 0.005)
    x0 = complex([0, 100, 2e4]); # initial values
    (d, tspan) = (1.0e-16, (0.0,1000)); # step size and time interval in days
    solution = sensitivity(x0, p, d, tspan); # find solution and partials

    (d, tspan) = (1.0e-16, (0.0,100)); # step size and time interval in days
    p1 = complex([0.0005, 0.1]); # kon koff (afinidad = 0.005)
    solution1 = sensitivity(x0, p1, d, tspan); # find solution and partials

    (d, tspan) = (1.0e-16, (0.0,10000)); # step size and time interval in days
    p2 = complex([5e-6, 0.001]); # kon koff (afinidad = 0.005)
    solution2 = sensitivity(x0, p2, d, tspan); # find solution and partials   

    ss1 = log10.(abs.(solution[1][:, 3]))
    ss2 = log10.(abs.(solution1[1][:, 3]))
    ss3 = log10.(abs.(solution2[1][:, 3]))
    pTlog = plot(ss1, label = "koff = 0.01", xlabel= "t", ylabel = "R(t)")
    plot!(ss2, label = "koff = 0.1")
    plot!(ss3, label = "koff = 0.001")
    display(pTlog)

    pT = plot(solution[1][:, 3], label = "koff = 0.01", xlabel= "t", ylabel = "R(t)")
    plot!(solution1[1][:, 3], label = "koff = 0.1")
    plot!(solution2[1][:, 3], label = "koff = 0.001")
    #Plots.plot(solution[1][:, 3], label = "x1", xlabel= "t", ylabel = "R(t)") #xlims = (tspan[1],tspan[2]))
    display(pT)

    p1 = plot(solution[1][:, 3], label = "koff = 0.01", xlabel= "t", ylabel = "R(t)", title = "Solución 1")
    display(p1)
    p2 = plot(solution1[1][:, 3], label = "koff = 0.1", xlabel= "t", ylabel = "R(t)", title = "Solución 2")
    display(p2)
    p3 = plot(solution2[1][:, 3], label = "koff = 0.001", xlabel= "t", ylabel = "R(t)", title = "Solución 3")
    display(p3)


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

    p = complex([5e-5, 0.01, 1]); # kon koff kp
    x0 = complex([100, 2e4, 0, 0]); # initial values
    (d, tspan) = (1.0e-16, (0.0,50)); # step size and time interval in days
    solution = sensitivity(x0, p, d, tspan); # find solution and partials
    Plots.plot(solution[4][:, 3], label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))

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
    Plots.plot(solution[13][:, 3], label = "x1", xlabel= "t", ylabel = "S") #xlims = (tspan[1],tspan[2]))

end

