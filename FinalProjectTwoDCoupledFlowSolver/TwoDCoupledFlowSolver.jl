#=
TwoDCoupledFlowSolver.jl
Jacob Child
April 11th, 2024
High Level Psuedocode: implement the SIMPLE CFD algorithm and solve for 2D channel flow
Problem: A 2D channel flow is to be solved using the SIMPLE CFD algorithm. Use an upwind scheme backward staggered grid with five pressure nodes and four velocity nodes. The stagnation pressur eis given at the inlet and the static pressure is specified at the exit. 
Assumptions: steady and viscous and the density of the fluid is constant.

Steps To Solve: Discritize, initialize, then calculate the new u using the A and b matrix maker and solver, then do the same for the pressure. 

=#
using Plots
include("TwoDSimpleFuncs.jl")

# Givens
ρ = 1000 #kg/m^3
μ = 0.001 #Pa*s
L = 5 / 100 #m 
W = 1 / 100 #m
#Boundary Conditions, no slip walls on the top and bottom
Vin = 0.001 #m/s
Pe = 0 #Pa
mdotin = ρ * Vin * W
#Initial guesses
u = 0.001 #m/s
v = 0.0001 #m/s
P  = 0.001 #Pa
#Under Relaxation Factors
αu = .5 #under relaxation factor for velocity, .4 kind of works for 64
αv = .5 #under relaxation factor for velocity
αp = 1 #under relaxation factor for pressure

#Discritize, 4x4 interior Pressure nodes, 5x4 interior u, 4x3 interior v
N = 64
P = ones(N,N) * P 
u = ones(N,N+1) * u
v = ones(N+1,N+2) * v
dx = L/N
dy = W/N

#=
#just iterate for a few times in a for loop 
us, vs, Ps = [], [], []
push!(us, u)
push!(vs, v)
push!(Ps, P)
for i in 1:10
    u3, v3, P3 = iterator(us[end],vs[end],Ps[end])
    push!(us, u3)
    push!(vs, v3)
    push!(Ps, P3)
end
=#
#us, vs, Ps = Converger(u,v,P)

#plot prep 
Pgx = range(dx/2,L-dx/2,step=dx)
Pgy = range(dy/2,W-dy/2,step=dy)
ugx = range(0, L, step=dx)
ugy = range(dy/2, W-dy/2, step=dy)
vgx = range(0-dx/2, L+dx/2, step=dx)
vgy = range(0, W, step=dy)
#for plotting
Pxcoords = repeat(Pgx, length(Pgy))
Pycoords = repeat(Pgy, inner = length(Pgx))
ucoordsx = repeat(ugx, length(ugy))
ucoordsy = repeat(ugy, inner = length(ugx))
vcoordsx = repeat(vgx, length(vgy))
vcoordsy = repeat(vgy, inner = length(vgx))
channelx = [0,0,L,L,0]
channely = [0,W,W,0,0]
GeoPlot = plot(channelx,channely,label = "channel", title="Geometry Discritization Plot")
scatter!(Pxcoords, Pycoords, label = "Pressure Nodes", markershape=:star5)
scatter!(ucoordsx, ucoordsy, label = "u Nodes", color =:red)
scatter!(vcoordsx, vcoordsy, label = "v Nodes",markershape =:cross, color =:black)
#plot the results, put us, vs, and Ps on individual contour plots
#final u plot in a filled contour plot including the channel
uPlot = contourf(ugx, ugy, us[end], title="u velocity Contour Plot (N = $N)", color=:viridis, xlabel="x (m)", ylabel="y (m)", aspect_ratio=:equal)
# scatter!(Pxcoords, Pycoords, label = "Pressure Nodes", markershape=:star5)
# scatter!(ucoordsx, ucoordsy, label = "u Nodes", color =:red)
# scatter!(vcoordsx, vcoordsy, label = "v Nodes",markershape =:cross, color =:black)
#final v plot in a filled contour plot including the channel
vPlot = contourf(vgx, vgy, vs[end], title="v velocity Contour Plot(N = $N)", color=:viridis, xlabel="x (m)", ylabel="y (m)", aspect_ratio=:equal)
#final P plot in a filled contour plot including the channel
PPlot = contourf(Pgx, Pgy, Ps[end], title="Pressure Contour Plot (N = $N)", color=:viridis, xlabel="x (m)", ylabel="y (m)", aspect_ratio=:equal)


#plot of pressure vs x along the centerline
midpoints = Int64.([N/2, N/2 + 1]) #assuming always an even N
Pcl = sum(Ps[end][midpoints,:],dims=1) / 2
PclPlot = plot(Pgx, vec(Pcl), title="Pressure vs x along the centerline (N = $N)", xlabel="x (m)", ylabel="Pa", label = "Pressure");

savefig(GeoPlot,"FinalProjectTwoDCoupledFlowSolver/GeoPlot.png")
savefig(uPlot,"FinalProjectTwoDCoupledFlowSolver/uPlot.png")
savefig(vPlot,"FinalProjectTwoDCoupledFlowSolver/vPlot.png")
savefig(PPlot,"FinalProjectTwoDCoupledFlowSolver/PPlot.png")
savefig(PclPlot,"FinalProjectTwoDCoupledFlowSolver/PclPlot.png")
writedlm("FinalProjectTwoDCoupledFlowSolver/Final64Pressure.csv", Ps[end], ',')
writedlm("FinalProjectTwoDCoupledFlowSolver/Final64U.csv", us[end], ',')
writedlm("FinalProjectTwoDCoupledFlowSolver/Final64V.csv", vs[end], ',')
writedlm("FinalProjectTwoDCoupledFlowSolver/Final64Pcl.csv", Pcl, ',')

PclCompPlot = plot(Pgx, vec(Pcl), title="Pressure vs x along the centerline", xlabel="x (m)", ylabel="Pa", label = "Coded 2D Flow Solver Pressure (N=64)");
plot!(CenterlinePressure01cfd[:,1],CenterlinePressure01cfd[:,2], label = "CFD Pressure")
savefig(PclCompPlot,"FinalProjectTwoDCoupledFlowSolver/PclCompPlot.png")

#plot of u velocity at the midline vs y 
midline = Int64(N/2)
uMidline = us[end][:, midline]
uMidlineCompPlot = plot(ugy, uMidline, title="u Velocity at the midline (x = 0.025m)", xlabel="y (m)", ylabel="m/s", label = "Coded 2D Flow Solver u Velocity (N=64)");
plot!(uVlmidline01cfd[:,1],uVlmidline01cfd[:,2], label = "CFD u Velocity", legend=:bottomright)
savefig(uMidlineCompPlot,"FinalProjectTwoDCoupledFlowSolver/uMidlineCompPlot.png")

#plot of v velocity at 0.0015 above the bottom wall vs x
ProfileLoc = argmin(abs.(vgy .- 0.0015))
vProfile64 = vs[end][end-ProfileLoc,:]
vProfileCompPlot = plot(vgx, vProfile64, title="v Velocity at 0.0015 above the bottom wall", xlabel="x (m)", ylabel="m/s", label = "Coded 2D Flow Solver v Velocity (N=64)");
plot!(Vlp01[:,1],Vlp01[:,2], label = "CFD v Velocity")
savefig(vProfileCompPlot,"FinalProjectTwoDCoupledFlowSolver/vProfileCompPlot.png")


#!see the next line
throw("stop here") 



#Grid independence study. Compare pressure vs x along centerline, and v velocity profile 
ProfileLoc = argmin(abs.(vgy .- 0.0015)) #ie whatever node is closest to 0.001
#initialize 
Pcls = []
vProfiles = []
uProfiles = []
for N in [4, 8, 16, 32]
    u = 0.001 #m/s
    v = 0.0001 #m/s
    P  = 0.001 #Pa
    #Under Relaxation Factors
    global αu = .5 #under relaxation factor for velocity
    global αv = .5 #under relaxation factor for velocity
    global αp = 1 #under relaxation factor for pressure
    dx = L/N
    dy = W/N
    P = ones(N,N) * P 
    u = ones(N,N+1) * u
    v = ones(N+1,N+2) * v
    us1, vs1, Ps1 = Converger(u,v,P)
    midpoints = Int64.([N/2, N/2 + 1]) #assuming always an even N
    push!(Pcls, sum(Ps1[end][midpoints,:],dims=1) / 2)
    midline1 = Int64(N/2)
    push!(uProfiles, us1[end][:, midline1])
    vgy1 = range(0, W, step=dy)
    ProfileLoc1 = argmin(abs.(vgy1 .- 0.0015))
    push!(vProfiles, vs1[end][end-ProfileLoc1,:])
end
# push!(Pcls, Pcl)
# push!(uProfiles, uMidline)
# push!(vProfiles, vProfile64)
PclConvPlot = plot(title="Convergence study of centerline pressure", xlabel="m", ylabel="Pa")
#for loop isn't working, so I will just do it manually
dx = L/4
Pgx = range(dx/2,L-dx/2,step=dx)
plot!(Pgx, vec(Pcls[1]), label = "N = 4")
dx = L/8
Pgx = range(dx/2,L-dx/2,step=dx)
plot!(Pgx, vec(Pcls[4]), label = "N = 8")
dx = L/16
Pgx = range(dx/2,L-dx/2,step=dx)
plot!(Pgx, vec(Pcls[5]), label = "N = 16")
dx = L/32
Pgx = range(dx/2,L-dx/2,step=dx)
plot!(Pgx, vec(Pcls[6]), label = "N = 32")
dx = L/64
Pgx = range(dx/2,L-dx/2,step=dx)
plot!(Pgx, vec(Pcls[7]), label = "N = 64")
plot!(title = "Convergence Study of Centerline Pressure", xlabel="x (m)", ylabel="Pa")
savefig(PclConvPlot,"FinalProjectTwoDCoupledFlowSolver/PclConvPlot.png")

vProfileConvPlot = plot(title="convergence study of v velocity Profile", xlabel="m/s", ylabel="m/s")
#for loop isn't working, so I will just do it manually
dx = L/4
vgx = range(0-dx/2, L+dx/2, step=dx)
plot!(vgx, vec(vProfiles[1]), label = "N = 4")
dx = L/8
vgx = range(0-dx/2, L+dx/2, step=dx)
plot!(vgx, vec(vProfiles[2]), label = "N = 8")
dx = L/16
vgx = range(0-dx/2, L+dx/2, step=dx)
plot!(vgx, vec(vProfiles[3]), label = "N = 16")
dx = L/32   
vgx = range(0-dx/2, L+dx/2, step=dx)
plot!(vgx, vec(vProfiles[4]), label = "N = 32")
dx = L/64
vgx = range(0-dx/2, L+dx/2, step=dx)
plot!(vgx, vec(vProfiles[5]), label = "N = 64")
plot!(title = "Convergence Study of v Velocity Profile", xlabel="x (m)", ylabel="m/s")
savefig(vProfileConvPlot,"FinalProjectTwoDCoupledFlowSolver/vProfileConvPlot.png")

uProfileConvPlot = plot(title="Convergence Study of u Velocity Profile", xlabel="y (m)", ylabel="m/s")
#for loop isn't working, so I will just do it manually
dy = W/4
ugy = range(dy/2, W-dy/2, step=dy)
plot!(ugy, vec(uProfiles[1]), label = "N = 4")
dy = W/8
ugy = range(dy/2, W-dy/2, step=dy)  
plot!(ugy, vec(uProfiles[2]), label = "N = 8")
dy = W/16
ugy = range(dy/2, W-dy/2, step=dy)
plot!(ugy, vec(uProfiles[3]), label = "N = 16")
dy = W/32
ugy = range(dy/2, W-dy/2, step=dy)
plot!(ugy, vec(uProfiles[4]), label = "N = 32")
dy = W/64
ugy = range(dy/2, W-dy/2, step=dy)
plot!(ugy, vec(uProfiles[5]), label = "N = 64")
plot!(title = "Convergence Study of u Velocity Profile", xlabel="y (m)", ylabel="m/s")
savefig(uProfileConvPlot,"FinalProjectTwoDCoupledFlowSolver/uProfileConvPlot.png")

for (i, (n, pcl)) in enumerate(zip([4, 8, 16, 32, 64], Pcls))
    dx = L/n
    Pgx = range(dx/2,L-dx/2,step=dx)
    plot!(Pgx, vec(pcl), label = "N = $(2^(i+1))")
end
#find the maximum value in each row of Pcl 
maxPcl = [maximum(row) for row in Pcls]
PclPercentdiff = [abs.(maxPcl[i] - maxPcl[i+1]) / maxPcl[i+1] * 100 for i in 1:length(maxPcl)-1]

#find the maximum value in each row of uProfiles
maxuProfiles = [maximum(row) for row in uProfiles]
uProfilePercentdiff = [abs.(maxuProfiles[i] - maxuProfiles[i+1]) / maxuProfiles[i+1] * 100 for i in 1:length(maxuProfiles)-1]

#find the maximum value in each row of vProfiles
maxvProfiles = [maximum(row) for row in vProfiles]
vProfilePercentdiff = [abs.(maxvProfiles[i] - maxvProfiles[i+1]) / maxvProfiles[i+1] * 100 for i in 1:length(maxvProfiles)-1]





#=
#Below is all of the code used for development, it has now all been put in functions.
#=

 #no need for plotting anymore
Pgx = range(dx/2,L-dx/2,step=dx)
Pgy = range(dy/2,W-dy/2,step=dy)
ugx = range(0, L, step=dx)
ugy = range(dy/2, W-dy/2, step=dy)
vgx = range(0-dx/2, L+dx/2, step=dx)
vgy = range(0, W, step=dy)
#for plotting
Pxcoords = repeat(Pgx, length(Pgy))
Pycoords = repeat(Pgy, inner = length(Pgx))
ucoordsx = repeat(ugx, length(ugy))
ucoordsy = repeat(ugy, inner = length(ugx))
vcoordsx = repeat(vgx, length(vgy))
vcoordsy = repeat(vgy, inner = length(vgx))
GeoPlot = plot([0,0,L,L,0],[0,W,W,0,0],label = "channel", title="Geometry Discritization Plot")
scatter!(Pxcoords, Pycoords, label = "Pressure Nodes", markershape=:star5)
scatter!(ucoordsx, ucoordsy, label = "u Nodes", color =:red)
scatter!(vcoordsx, vcoordsy, label = "v Nodes",markershape =:cross, color =:black);
=#

#Begin to solve, to start I will do 1 iteration at a time
#preallocate us
aues = zero(u) #au easts
auws = zero(u) #au wests
auns = zero(u)
auss = zero(u) 
aups = zero(u)
Bus = zero(u)
#preallocate vs 
aves = zero(v) #av easts
avws = zero(v) #av wests
avns = zero(v)
avss = zero(v)
avps = zero(v)
Bvs = zero(v)
#Preallocate Ps 
apes = zero(P) #ap easts
apws = zero(P) #ap wests
apns = zero(P)
apss = zero(P)
apps = zero(P)
Bps = zero(P)


#upper wall 
aues[1,:], auws[1,:], auns[1,:], auss[1,:], aups[1,:], Bus[1,:] = uABMaker(u[1:3,:], v[1:2,:], P[1,:], wallflag = "upper")

#internal  for loop only 
for i in axes(u[1:end-2,:],1)
aues[i+1,:], auws[i+1,:], auns[i+1,:], auss[i+1,:], aups[i+1,:], Bus[i+1,:]  = uABMaker(u[i:i+2,:], v[i+1:i+2,:], P[i+1,:])
end

#bottom wall
aues[end,:], auws[end,:], auns[end,:], auss[end,:], aups[end,:], Bus[end,:] = uABMaker(u[end-2:end,:], v[end-1:end,:], P[end,:], wallflag = "lower")

Au = AMaker2D(reverse(aups,dims=1),reverse(aues,dims=1),reverse(auws,dims=1),reverse(auns,dims=1),reverse(auss,dims=1),size(u))
u1 = reverse(reshape(Au \ vec(reverse(Bus,dims=1)'),reverse(size(u)))',dims=1)
#mdot correction
mdotinapparent = sum(u1[:,1]) 
mdotoutappar = sum(u1[:,end])
u1[:,end] = u1[:,end] * mdotinapparent / mdotoutappar

#v calculations
#upper wall
aves[1,:], avws[1,:], avns[1,:], avss[1,:], avps[1,:], Bvs[1,:] = vABMaker(u1[1:2,:], v[1:3,:], P[1:2,:], wallflag = "upper")
#internal  for loop only 
for i in axes(v[1:end-2,:],1)   
    aves[i+1,:], avws[i+1,:], avns[i+1,:], avss[i+1,:], avps[i+1,:], Bvs[i+1,:] = vABMaker(u1[i:i+1,:], v[i:i+2,:], P[i:i+1,:])

end

#bottom wall
aves[end,:], avws[end,:], avns[end,:], avss[end,:], avps[end,:], Bvs[end,:] = vABMaker(u1[end-1:end,:], v[end-2:end,:], P[end-1:end,:], wallflag = "lower")

Av = AMaker2D(reverse(avps,dims=1),reverse(aves,dims=1),reverse(avws,dims=1),reverse(avns,dims=1),reverse(avss,dims=1),size(v))
v1 = reverse(reshape(Av \ vec(reverse(Bvs,dims=1)'),reverse(size(v)))', dims=1)

#P calculations
#upper wall
apes[1,:], apws[1,:], apns[1,:], apss[1,:], apps[1,:], Bps[1,:] = pPrimeABMaker(aups[1,:],avps[1:2,:],u1[1,:],v1[1:2,:],P[1,:], wallflag = "upper")
#internal  for loop only
for i in axes(P[1:end-2,:],1)
apes[i+1,:], apws[i+1,:], apns[i+1,:], apss[i+1,:], apps[i+1,:], Bps[i+1,:] = pPrimeABMaker(aups[i+1,:],avps[i+1:i+2,:],u1[i+1,:],v1[i+1:i+2,:],P[i+1,:])
end
#bottom wall 
apes[end,:], apws[end,:], apns[end,:], apss[end,:], apps[end,:], Bps[end,:] = pPrimeABMaker(aups[end,:],avps[end-1:end,:],u1[end,:],v1[end-1:end,:],P[end,:], wallflag = "lower")

AP = AMaker2D(reverse(apps,dims=1),reverse(apes,dims=1),reverse(apws,dims=1),reverse(apns,dims=1),reverse(apss,dims=1),size(P))
P1 = reverse(reshape(AP \ vec(reverse(Bps,dims=1)'),reverse(size(P)))', dims=1)


#Calculate the new P, u, v
Pnew = P .+ αp .* P1
diJ = dy ./ aups
dIj = dx ./ avps
unew = deepcopy(u1)
vnew = deepcopy(v1)
unew[:,2:end-1] = u1[:,2:end-1] .+ diJ[:,2:end-1] .* -diff(P1,dims=2)
vnew[2:4,2:5] = v1[2:4,2:5] .+ dIj[2:4,2:5] .* diff(P1,dims=1)
#enforce the boundary condition
vnew[:,end] .= 0

=#

