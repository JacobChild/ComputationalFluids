#=
HW4.jl
Feb 19th, 2024
Jacob Child
Pseudocode: Do HW3, but use Gauss-Seidel
=#

using Plots
include("HW4Funcs.jl")

#Problem 1
L = .02
D = .003
ϵ = 0
σ = 5.67*10^-8
K = 401
h = 10
Tb = 400
Tinf = 273
Tsurr = 273
N = 80
Nodes, fw, fe, dx = PracticeB([0,L], N)
#Homogeneous material, so constant K and no qdot 
ks = ones(N+2) * K
Tguess = ones(N+2) * Tb
T = Tguess
#do a quick iterate 3 times
aw, ap, ae, b = abmaker(Tguess, ks, dx, D, ϵ, σ, h, Tb, Tinf, Tsurr)

#Gauss Seidel by hand first to understand
Test1 = GaussSeidel(Tguess, aw, ap, ae, b)
#Comparer = TDMASolver(ap, ae, aw, b)
#while loop to iterate until error is less than .0001
errors = []
Ts = []
push!(errors, 1)
push!(Ts, Test1)
while errors[end] > 1*10^-6
    Test2 = GaussSeidel(Ts[end], aw, ap, ae, b)
    push!(errors, maximum(abs.(Test2 - Ts[end])))
    push!(Ts, Test2)
end
println("The number of iterations to get the error below 1*10^-6 is ", length(errors))
#Redo the iterations, but now with heat flux
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[end-3] - Test1[end-2]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidel(Ts[end], aw, ap, ae, b)
    TipHeatFlux = -K * (Test2[end-2] - Test2[end-1]) / dx
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("The number of iterations to get the heat flux error below 1*10^-6 is ", length(errors))
#redo but with base heat flux 
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[2] - Test1[1]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidel(Ts[end], aw, ap, ae, b)
    TipHeatFlux = -K * (Test2[1] - Test2[2]) / (.5*dx)
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("The number of iterations to get the base heat flux error below 1*10^-6 is ", length(errors))

#Problem 2
ϵ = 1;
#N = 20
Nodes, fw, fe, dx = PracticeB([0,L], N)
#Homogeneous material, so constant K and no qdot 
ks = ones(N+2) * K
Tguess = ones(N+2) * Tb
T = Tguess
#do a quick iterate 3 times
aw, ap, ae, b = abmaker(Tguess, ks, dx, D, ϵ, σ, h, Tb, Tinf, Tsurr)
#Gauss Seidel by hand first to understand
Test1 = GaussSeidelRelaxed(Tguess, aw, ap, ae, b, 1)
#iterate based off of base heat flux convergence
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[2] - Test1[1]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidelRelaxed(Ts[end], aw, ap, ae, b, 1)
    TipHeatFlux = -K * (Test2[1] - Test2[2]) / (.5*dx)
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("Part 2 The number of iterations to get the base heat flux error below 1*10^-6 is ", length(errors))
N1 = length(errors)

#Gauss Seidel by hand first to understand
Test1 = GaussSeidelRelaxed(Tguess, aw, ap, ae, b, 1.2)
#iterate based off of base heat flux convergence
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[2] - Test1[1]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidelRelaxed(Ts[end], aw, ap, ae, b, 1.2)
    TipHeatFlux = -K * (Test2[1] - Test2[2]) / (.5*dx)
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("Part 2 The number of iterations with relaxation to get the base heat flux error below 1*10^-6 is ", length(errors))
N2 = length(errors)
percentdif = (N2 - N1) / N1 * 100
println("Part 2 The percent difference between the two is ", percentdif)

#Part 3
n = sqrt( 4*h / (D * K))
Tanalytical = @. cosh(n*(L - Nodes)) / cosh(n*L) * (Tb - Tinf) + Tinf

#Gauss Seidel by hand first to understand
Test1 = Tanalytical
#iterate based off of base heat flux convergence
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[2] - Test1[1]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidelRelaxed(Ts[end], aw, ap, ae, b, 1)
    TipHeatFlux = -K * (Test2[1] - Test2[2]) / (.5*dx)
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("Part 3 The number of iterations to get the base heat flux error below 1*10^-6 is ", length(errors))
N1 = length(errors)

#Gauss Seidel by hand first to understand
Test1 = Tanalytical
#iterate based off of base heat flux convergence
errors = []
qs = []
Ts = []
push!(Ts, Test1)
push!(errors, 1)
push!(qs, -K * (Test1[2] - Test1[1]) / dx)
while errors[end] > 1*10^-6
    Test2 = GaussSeidelRelaxed(Ts[end], aw, ap, ae, b, 1.2)
    TipHeatFlux = -K * (Test2[1] - Test2[2]) / (.5*dx)
    push!(errors, abs(qs[end] - TipHeatFlux))
    push!(Ts, Test2)
    push!(qs, TipHeatFlux)
end
println("Part 3 The number of iterations with relaxation to get the base heat flux error below 1*10^-6 is ", length(errors))
N2 = length(errors)
percentdif = (N2 - N1) / N1 * 100
println("Part 3 The percent difference between the two is ", percentdif)
