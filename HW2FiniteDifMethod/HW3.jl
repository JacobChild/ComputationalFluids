#= 
HW3.jl
Feb 13th, 2024
Jacob Child
=#
using Plots


#Problem 1
#Same as before, but Kb is a function of x, so I will use MaterialMaker5 and EasyRun5 to solve it.
function problem1()

    bounds = [0,.05]
    #Material properties
    ka = 15 #W/(m*K)
    kb = 60 #W/(m*K)
    qdota = 4*10^6 #W/m^3
    h = 1000 #W/(m^2*K)
    T∞ = 300 #K
    MatInterfaceLoc = .03 #m
    Kbfoo(xf) = 137*exp(25*xf-2) #Kb as a function of x
    N = 5
    T5, _, _ = EasyRun5(bounds, N, ka, Kbfoo, qdota, h, T∞, MatInterfaceLoc)
    N = 10
    T10, _, _ = EasyRun5(bounds, N, ka, Kbfoo, qdota, h, T∞, MatInterfaceLoc)
    N = 20
    T20, _, _ = EasyRun5(bounds, N, ka, Kbfoo, qdota, h, T∞, MatInterfaceLoc)
    #Problem 1 Richardson Extrapolation
    #Part A 
    #Estimate the order of the simulation used 
    P = log((T10[1] - T5[1])/ (T20[1] - T10[1])) / log(2)
    println("The order of the simulation is $(round(P,digits=2)).")
    Ts = [T5[1], T10[1], T20[1]]
    println("For Ts = $Ts")
    #Part B 
    #Estimate the Grid Converged value
    TFinal = T20[1] + (T10[1] - T20[1]) / (1 - 2^(round(P,digits=2)))
    println("The grid converged value is $(round(TFinal,digits=2)).")
    #Part C
    #Plot the base temperature T[1] for 5,10,20,40,5120,10240 CVs and show the TFinal from part b on this plot 
    N = [5,10,20,40,5120,10240]
    #PartCT = zero(Float64.(N))
    #=
    for (i,cn) in enumerate(N)
        TCurrent, _, _ = EasyRun5(bounds, cn, ka, Kbfoo, qdota, h, T∞, MatInterfaceLoc)
    end
    =#
    Prb1Plt = plot(N, PartCT, xscale=:log10, label = "T[1] for N CVs", xlabel = "Number of Control Volumes", ylabel = "Temperature (K)", title = "Boundary Temp vs # of Control Volumes", marker=:x)
    hline!([TFinal], label = "Predicted TFinal from Richardson Extrapolation", style = :dash, color = :red)
    #Generate a Table that lists the number of control volumes, the temperature at the left boundary, and the percent difference between the temperature at the left boundary and the predicted TFinal
    PercentDiff = (PartCT .- TFinal) ./ TFinal * 100
    Table = [N PartCT round.(PercentDiff,digits=4)]
    TableHeader = ["N" "T[1]" "Percent Difference"]
    FullTable = [TableHeader; Table]
    return Prb1Plt, FullTable
end

#Problem 2
L = .02
D = .003
ϵ = 0
σ = 5.67*10^-8
K = 401
h = 10
Tb = 400
Tinf = 273
Tsurr = 273
N = 50
Nodes, fw, fe, dx = PracticeB([0,L], N)
#Homogeneous material, so constant K and no qdot 
ks = ones(N+2) * K
Tguess = ones(N+2) * Tb
T = Tguess .- 1
#do a quick iterate 3 times
aw, ap, ae, b = abmaker(Tguess, ks, dx, D, ϵ, σ, h, Tb, Tinf, Tsurr)
T1 = TDMASolver(ap, ae, aw, b)
Tguess = T1
#analytical solution
n = sqrt( 4*h / (D * K))
Tanalytical = @. cosh(n*(L - Nodes)) / cosh(n*L) * (Tb - Tinf) + Tinf
#iterate to find needed CVs to have max error less than .0001 by comparing against Tanalytical
Ts = []
Ns = []
errors = []
push!(errors, 1)
push!(Ts, T1)
push!(Ns, N)
#throw("whoa tiger, not ready yet")
while  errors[end] > .0001
    Nodesw, _, _, dxw = PracticeB([0,L], Ns[end])
    ksw = ones(Ns[end]+2) * K
    Tguessw = ones(Ns[end]+2) * Tb
    aww, apw, aew, bw = abmaker(Tguessw, ksw, dxw, D, ϵ, σ, h, Tb, Tinf, Tsurr)
    Tw = TDMASolver(apw, aew, aww, bw)
    Tanalyticalw = @. cosh(n*(L - Nodesw)) / cosh(n*L) * (Tb - Tinf) + Tinf
    # Check the lengths of Tw and Tanalyticalw
    push!(errors, maximum(abs.(Tw - Tanalyticalw)))
    push!(Ts, Tw)
    push!(Ns, Ns[end].+1)
end
#plot
Prb2Plt = plot(Ns, errors, yscale=:log10, label = "Error", xlabel = "Number of Control Volumes", ylabel = "Error", title = "Error vs # of Control Volumes", marker=:x)
Ttip2 = Ts[end-1][end] #this is with 50 CVs


#Problem 3
ϵ = 1
h = 10
Ts3 = []
errors = []
#see if I can mess it up, 
Tguess = ones(N+2)
push!(Ts3, Tguess)
push!(errors, 1)
#iterate, keeping N the same, but changing the input T 

while  errors[end] > .0001
    aww, apw, aew, bw = abmaker(Ts3[end], ks, dx, D, ϵ, σ, h, Tb, Tinf, Tsurr)
    Tgw = TDMASolver(apw, aew, aww, bw)
    push!(errors, maximum(abs.(Tgw - Ts3[end])))
    push!(Ts3, Tgw)
end
Ttip3 = Ts3[end][end] 

Prb3Plt = plot(Ts3[end], label = "Final T", xlabel = "Node", ylabel = "Temperature (K)", title = "Final T vs Node", marker=:x)
#for fun, overlay a heat map of the temperature
# Assuming Ts3[end] is your 1D temperature distribution
temperature = Ts3[end]
# Create a new array with the same size as temperature, filled with 1's
# This is used to create a 2D plot from the 1D temperature data
y = ones(length(temperature)) * 399.25
# Create a plot with color mapped to the temperature
plot!(y, linewidth = 20, line_z = temperature, label = "Rod", c=:plasma)