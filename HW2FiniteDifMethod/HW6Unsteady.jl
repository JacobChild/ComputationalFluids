#=
HW6Unsteady.jl
March 19th, 2024
Jacob Child
=#
using Plots, LaTeXStrings
include("HW6Funcs.jl")

#Givens
w = 1
L = 10 / 1000 #m (10mm)
bthick = 0.3 / 1000 #m (.3mm)
ρ = 2770 #kg/m^3
c = 875 #Jacob
k = 177 #W/mK
T∞ = 293 #K
Tsurr = 293 #K
Tinit = 293 #K
ϵ = [0,1]
h = [5,30] #W/m^2K
σ = 5.67e-8 #W/m^2K^4
Tb = 353 #K

#Discretization (Practice A)
N = 100
Nodes, CS, dx = PracticeA([0,L], N) #remember dx is a vector of length N with the different dx values, so I do not! need to adjust in the a,b maker
#Solve 
Tsh5e0 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[1], ϵ[1], σ, Tb, ρ, c, k)
Tsh5e1 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[1], ϵ[2], σ, Tb, ρ, c, k)
Tsh30e0 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[2], ϵ[1], σ, Tb, ρ, c, k)
Tsh30e1 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[2], ϵ[2], σ, Tb, ρ, c, k)

#Perform a grid convergence study
N = [25, 50, 100, 200, 400, 800]
# I only want the time portion of the loop to run once, so 
tmax = 0.1 #s
Δt = .1 #s
#for loop to iterate through N, compare only the h=30, ϵ=1 case as it should be the most dramatically changing, save and compare base heat flux as shown in HW4
Qbase = zeros(length(N))
for (i, Nf) in enumerate(N)
    Nodes, CS, dx = PracticeA([0,L], Nf)
    T = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[2], ϵ[2], σ, Tb, ρ, c, k)
    Qbase[i] = -k*(T[end][2] - T[end][1])/dx[2]
end
#Richardson Extrapolation
#estimate the order of the simulation used 
P = log((Qbase[3] - Qbase[2])/ (Qbase[4] - Qbase[3])) / log(2)
println("The order of the simulation is $(round(P,digits=2)).")
#Estimate the Grid Converged value
QFinal = Qbase[4] + (Qbase[3] - Qbase[4]) / (1 - 2^(round(P,digits=2)))
println("The grid converged value is $(round(QFinal,digits=2)).")
GridConvergencePlot = plot(N, Qbase, label = "Q at the base", xlabel = "Number of Control Volumes", ylabel = L"Q (\frac{W}{m^2})", title = "Q at the base vs # of Control Volumes", marker=:x)
hline!([QFinal], label = "Predicted QFinal from Richardson Extrapolation", style = :dash, color = :red)
#Generate a Table that lists the number of control volumes and the percent difference between the heat flux at the base and the predicted QFinal
PercentDiff = (Qbase .- QFinal) ./ QFinal * 100
Table = [N round.(PercentDiff,digits=2)]
TableHeader = ["N" "Percent Difference"]
GridConvergenceTable = [TableHeader; Table]

#Perform a timestep indpendence study
N = 200
Δt = [.1, .05, .025, .0125, .00625, .003125]
tmax = 5 #s
#for loop to iterate through Δt, compare only the h=30, ϵ=1 case as it should be the most dramatically changing, save and compare base heat flux as shown in HW4
Qbaset = zeros(length(Δt))
for (i, Δtf) in enumerate(Δt)
    Nodes, CS, dx = PracticeA([0,L], N)
    T = TimeStepper(Nodes, dx, Tinit, Δtf, tmax, T∞, Tsurr, bthick, h[2], ϵ[2], σ, Tb, ρ, c, k)
    Qbaset[i] = -k*(T[end][2] - T[end][1])/dx[2]
end
#? I'm not sure you can do Richardson Extrapolation with time stepping. The results did not seem to converge to the RE predicted value, so I am commenting it out.
#estimate the order of the simulation used
#P = log((Qbaset[3] - Qbaset[2])/ (Qbaset[4] - Qbaset[3])) / log(2)
#println("The order of the simulation is $(round(P,digits=2)).")
#Estimate the Grid Converged value
#QFinalt = Qbaset[4] + (Qbaset[3] - Qbaset[4]) / (1 - 2^(round(P,digits=2)))
#println("The timestep converged value is $(round(QFinalt,digits=2)).")
TimeStepConvergencePlot = plot(Δt, Qbaset, label = "Q at the base", xlabel = "Timestep (s)", ylabel = L"Q (\frac{W}{m^2})", title = "Q at the base vs Timestep", marker=:x, xflip=true, legend=:topright)
#hline!([QFinalt], label = "Predicted QFinal from Richardson Extrapolation", style = :dash, color = :red)
#Generate a Table that lists the timestep and the percent difference between the heat flux at the base at each timestep and the one at the next timestep
PercentDifft = (Qbaset[1:end-1] .- Qbaset[2:end]) ./ Qbaset[2:end] .* 100
pushfirst!(PercentDifft, 0)
Tablet = [Δt round.(PercentDifft,digits=2)]
TableHeadert = ["Δt" "Percent Change"]
TimeStepConvergenceTable = [TableHeadert; Tablet]


#Plot Q(t) at the base for the 4 cases
#Use the good values after the convergence study 
Δt = .00625 #s
tmax = 5 #s
N = 200 
times = range(start=0, stop=tmax, step=Δt)
Nodes, CS, dx = PracticeA([0,L], N) #remember dx is a vector of length N with the different dx values, so I do not! need to adjust in the a,b maker
#Solve 
Tsh5e0 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[1], ϵ[1], σ, Tb, ρ, c, k)
Tsh5e1 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[1], ϵ[2], σ, Tb, ρ, c, k)
Tsh30e0 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[2], ϵ[1], σ, Tb, ρ, c, k)
Tsh30e1 = TimeStepper(Nodes, dx, Tinit, Δt, tmax, T∞, Tsurr, bthick, h[2], ϵ[2], σ, Tb, ρ, c, k)
#comprehension through Ts at the base (first two nodes) to get Q at the base through time 
Qsh5e0 = [-k*(Tsh5e0[i][2] - Tsh5e0[i][1])/dx[2] for i in 1:length(Tsh5e0)]
Qsh5e1 = [-k*(Tsh5e1[i][2] - Tsh5e1[i][1])/dx[2] for i in 1:length(Tsh5e1)]
Qsh30e0 = [-k*(Tsh30e0[i][2] - Tsh30e0[i][1])/dx[2] for i in 1:length(Tsh30e0)]
Qsh30e1 = [-k*(Tsh30e1[i][2] - Tsh30e1[i][1])/dx[2] for i in 1:length(Tsh30e1)]
#plot 
Qplot = plot(times, Qsh5e0, label = "h = 5, ϵ = 0", xlabel = "Time (s)", ylabel = L"Q (\frac{W}{m^2}, log scale)", title = "Q at the base through time", marker=:x, markeralpha = 0.25)
plot!(times, Qsh5e1, label = "h = 5, ϵ = 1", marker=:circle, markeralpha = 0.25)
plot!(times, Qsh30e0, label = "h = 30, ϵ = 0", marker=:square, markeralpha = 0.25)
plot!(times, Qsh30e1, label = "h = 30, ϵ = 1", marker=:diamond, markeralpha = 0.25)
plot!(yscale=:log)

#Part 4 compare with the exact solution
x = range(0, stop=L, length=N) # Define the x values (from 0 to L)
Ac = w*bthick
P = 2*w + 2*bthick
m = sqrt.(h.*P/(k*Ac))
Th5 = @. T∞ + (cosh(m[1]*(L-x)) / cosh(m[1]*L) ) * (Tb - T∞)
Th30 = @. T∞ + (cosh(m[2]*(L-x)) / cosh(m[2]*L) ) * (Tb - T∞)

#Plot the predicted profiles of the steady-state temperature distribution and compare to the exact solution 
TCompExactPlot = plot(x, Tsh5e0[end], label = "h = 5, ϵ = 0", xlabel = "x", ylabel = "Temperature", title = "Temperature distribution at t = 5s", marker=:x)
plot!(x, Tsh5e1[end], label = "h = 5, ϵ = 1", marker=:circle)
plot!(x, Tsh30e0[end], label = "h = 30, ϵ = 0", marker=:square)
plot!(x, Tsh30e1[end], label = "h = 30, ϵ = 1", marker=:diamond)
plot!(x, Th5, label = "Exact h = 5", linestyle = :dashdot, marker=:star)
plot!(x, Th30, label = "Exact h = 30", linestyle = :dashdot, marker=:star)


throw("Stop here")
# Create the animation
anim = @animate for i in 1:length(Ts)
    plot(x, Ts[i], title="Temperature distribution over time", xlabel="x", ylabel="Temperature", ylims = (Tinit, Tb))
end

# Save the animation as a GIF
gif(anim, "temperature.gif", fps = 10)
