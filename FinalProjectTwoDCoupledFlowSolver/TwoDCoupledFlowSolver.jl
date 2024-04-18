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
αu = .5 #under relaxation factor for velocity
αv = .5 #under relaxation factor for velocity
αp = 1 #under relaxation factor for pressure

#Discritize, 4x4 interior Pressure nodes, 5x4 interior u, 4x3 interior v
N = 4
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

