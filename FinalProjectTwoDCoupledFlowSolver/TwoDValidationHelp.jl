#=
TwoDValidationHelp.jl
Jacob Child
April 24th, 2024
Pseudocode: some quick calculations to help set up the star ccm validation of my 2D flow solver.
=#

using Plots, DelimitedFiles

h = 1/ 100 #m
l = 5/100 #m
ρ = 1000 #kg/m^3
μ = 0.001 #Pa*s
Vin = 0.001 #m/s
Pe = 0 #Pa

Re = ρ * Vin * l / μ #Re = 50, very low, safe to assume laminar

#Prism layer calculations
r = 1.2 #growth rate
H = 1.72*l/sqrt(Re) #blasius δstar for laminar flow, instead of schlicting for turbulent
#assume y+ = 1
cf = .664/sqrt(Re) #blasius for laminar flow, instead of schlicting for turbulent
h = l/Re *sqrt(2/cf) #height of the first cell
n = log(r,H/h*(r-1)+1) #number of cells, 14, (13.415)

#Convergence values
B = [.01,.02,.04,.08] #Base size
Pcl = [.005945, .005949, .005915, .005896] #Pressure at centerline entrance
Umc = [.001496, .001495, .001490, .001486] #u velocity at mid and centerline 
Vcle = [-.000010469, .000001125, .000009025, .0000066756] #v velocity at centerline exit
#Percent changes between values with the value before
PclPercChange = [(Pcl[i]-Pcl[i-1])/Pcl[i-1] * 100 for i in 2:length(Pcl)]
UmcPercChange = [(Umc[i]-Umc[i-1])/Umc[i-1] * 100 for i in 2:length(Umc)]
VclePercChange = [(Vcle[i]-Vcle[i-1])/Vcle[i-1] * 100 for i in 2:length(Vcle)]



#Plots
PercChangePlot = plot(B[1:end-1],PclPercChange, label="Pressure Centerline Percent Change", xlabel="Base Size", ylabel="Percent Change", title="Percent Change")
plot!(B[1:end-1],UmcPercChange, label="U Mid Centerline Percent Change")
#plot!(B[1:end-1],VclePercChange, label="V Centerline Exit Percent Change")

#compare v velocity at centerline profiles instead of specific value
Vcle01, header = readdlm("J:\\ComputationalFluids541\\CenterlineVVelocity01.csv", header = true, ',')
Vcle02, header = readdlm("J:\\ComputationalFluids541\\CenterlineVVelocity02.csv", header = true, ',')
Vcle04, header = readdlm("J:\\ComputationalFluids541\\CenterlineVVelocity04.csv", header = true, ',')
Vcle08, header = readdlm("J:\\ComputationalFluids541\\CenterlineVVelocity08.csv", header = true, ',')
VclePlot = plot(Vcle01[:,1],Vcle01[:,2], label="Base Size 0.01")
plot!(Vcle02[:,1],Vcle02[:,2], label="Base Size: 0.02m")
plot!(Vcle04[:,1],Vcle04[:,2], label="Base Size: 0.04m")
plot!(Vcle08[:,1],Vcle08[:,2], label="Base Size: 0.08m")
plot!(title = "Centerline v Velocity Profiles", xlabel = "x", ylabel = "v (m/s)")
#The above does not show any convergence, I am going to repeat, but not at the centerline, instead at 0.0015 above the bottom wall
#TODO: readdlm and plot the v velocity at 0.0015 above the bottom wall for the four base sizes
#TODO: export all of the u, v, P, contour plots to the to my doc
Vlp01, header = readdlm("J:\\ComputationalFluids541\\VVelocity0015NearEdge01.csv", header = true, ',')
Vlp02, header = readdlm("J:\\ComputationalFluids541\\VVelocity0015NearEdge02.csv", header = true, ',')
Vlp04, header = readdlm("J:\\ComputationalFluids541\\VVelocity0015NearEdge04.csv", header = true, ',')
Vlp08, header = readdlm("J:\\ComputationalFluids541\\VVelocity0015NearEdge08.csv", header = true, ',')
VlpPlot = plot(Vlp01[:,1],Vlp01[:,2], label="Base Size 0.01")
plot!(Vlp02[:,1],Vlp02[:,2], label="Base Size: 0.02m")
plot!(Vlp04[:,1],Vlp04[:,2], label="Base Size: 0.04m")
plot!(Vlp08[:,1],Vlp08[:,2], label="Base Size: 0.08m")
plot!(title = "v Velocity Profiles at 0.0015 above the bottom wall", xlabel = "x", ylabel = "v (m/s)")

#TODO: Pcl plot from 0.01 grid size to compare
uVlmidline01cfd, header = readdlm("J:\\ComputationalFluids541\\UVelocMidline01.csv", header = true, ',')
CenterlinePressure01cfd, header = readdlm("J:\\ComputationalFluids541\\CenterlinePressure01.csv", header = true, ',')