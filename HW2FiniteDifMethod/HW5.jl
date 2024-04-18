#=
HW5.jl
March 6th, 2024
Jacob Child
Pseudocode: Use practice A to discritize, then create the A matrix and b vector using one of the following schemes: Central, Upwind, Hybrid, Power Law. Compare the results of the different schemes in a table, then a plot.
=#
using Plots, XLSX, DataFrames
include("HW5Funcs.jl")

#Problem 1
bounds = [0.0,1.0]
N = 21
u = [1, 20, 75]
ρ = 1
Γ = 1
Nodes, CS, dx = PracticeA(bounds, N)
Errors = zeros(4,3)

P, F, D = PecletNum(Nodes, u[1], ρ, Γ, dx)
#central = 1, upwind = 2, hybrid = 3, power law = 4
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 1)
T1 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 2)
T2 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 3)
T3 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 4)
T4 = TDMASolver(Ap, Ae, Aw, b)
Panalytical = [ρ*u[1]*nodef/Γ for nodef in Nodes]
Pfull = @. ρ*u*bounds[2]/Γ
Tanalytical = @. 1 - (exp(Panalytical/bounds[2]) - 1)/(exp(Pfull[1]) - 1)
Headers = ["Central" "Upwind" "Hybrid" "Power Law" "Analytical"]
Table1 = [T1 T2 T3 T4 Tanalytical]
FullTable1 = [Headers; Table1]
XLSX.writetable("Table1.xlsx", DataFrame(FullTable1,:auto))
Errors[1,1] = sum(Tanalytical .- T1)
Errors[2,1] = sum(Tanalytical .- T2)
Errors[3,1] = sum(Tanalytical .- T3)
Errors[4,1] = sum(Tanalytical .- T4)
#Plot
U1Plot = plot(Nodes, T1, label = "Central", xlabel = "Node Location", ylabel = "Phi", title = "Phi vs Node Location at U = 1", marker=:x)
plot!(Nodes, T2, label = "Upwind", marker=:circle)
plot!(Nodes, T3, label = "Hybrid", marker=:square)
plot!(Nodes, T4, label = "Power Law", marker=:diamond)
plot!(Nodes, Tanalytical, label = "Analytical", linestyle = :dashdot, marker=:star);


#redo at next u 
P, F, D = PecletNum(Nodes, u[2], ρ, Γ, dx)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 1)
T1 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 2)
T2 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 3)
T3 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 4)
T4 = TDMASolver(Ap, Ae, Aw, b)
Panalytical = [ρ*u[2]*nodef/Γ for nodef in Nodes]
Pfull = @. ρ*u*bounds[2]/Γ
Tanalytical = @. 1 - (exp(Panalytical/bounds[2]) - 1)/(exp(Pfull[2]) - 1)
Table2 = [T1 T2 T3 T4 Tanalytical]
FullTable2 = [Headers; Table2]
XLSX.writetable("Table2.xlsx", DataFrame(FullTable2,:auto))
Errors[1,2] = sum(Tanalytical .- T1)
Errors[2,2] = sum(Tanalytical .- T2)
Errors[3,2] = sum(Tanalytical .- T3)
Errors[4,2] = sum(Tanalytical .- T4)
#Plot 
U2Plot = plot(Nodes, T1, label = "Central", xlabel = "Node Location", ylabel = "Phi", title = "Phi vs Node Location at U = 20", marker=:x)
plot!(Nodes, T2, label = "Upwind", marker=:circle)
plot!(Nodes, T3, label = "Hybrid", marker=:square)
plot!(Nodes, T4, label = "Power Law", marker=:diamond)
plot!(Nodes, Tanalytical, label = "Analytical", linestyle = :dashdot, marker=:star)

#redo at next u
P, F, D = PecletNum(Nodes, u[3], ρ, Γ, dx)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 1)
T1 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 2)
T2 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 3)
T3 = TDMASolver(Ap, Ae, Aw, b)
Ae, Aw, Ap, b = AbMakerHw5(P, F, D, dx, 4)
T4 = TDMASolver(Ap, Ae, Aw, b)
Panalytical = [ρ*u[3]*nodef/Γ for nodef in Nodes]
Pfull = @. ρ*u*bounds[2]/Γ
Tanalytical = @. 1 - (exp(Panalytical/bounds[2]) - 1)/(exp(Pfull[3]) - 1)
Table3 = [T1 T2 T3 T4 Tanalytical]
FullTable3 = [Headers; Table3]
XLSX.writetable("Table3.xlsx", DataFrame(FullTable3,:auto))
Errors[1,3] = sum(Tanalytical .- T1)
Errors[2,3] = sum(Tanalytical .- T2)
Errors[3,3] = sum(Tanalytical .- T3)
Errors[4,3] = sum(Tanalytical .- T4)
#Plot 
U3Plot = plot(Nodes, T1, label = "Central", xlabel = "Node Location", ylabel = "Phi", title = "Phi vs Node Location at U = 75", marker=:x)
plot!(Nodes, T2, label = "Upwind", marker=:circle)
plot!(Nodes, T3, label = "Hybrid", marker=:square)
plot!(Nodes, T4, label = "Power Law", marker=:diamond)
plot!(Nodes, Tanalytical, label = "Analytical", linestyle = :dashdot, marker=:star);

ErrorHeadings = ["Scheme" "1" "20" "75"]
ErrorSidebar = ["Central", "Upwind", "Hybrid", "Power Law"]
FullErrors = [ErrorHeadings; ErrorSidebar Errors]