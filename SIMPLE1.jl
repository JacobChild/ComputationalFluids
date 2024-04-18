#=
SIMPLE1.jl
Jacob Child
April 3rd, 2024
High Level Psuedocode: implement the SIMPLE CFD algorithm and solve for a planar 1D converging nozzle flow
Problem: A planar 1D converging nozzle flow is to be solved using the SIMPLE CFD algorithm. Use the backward staggered grid with five pressure nodes and four velocity nodes. The stagnation pressur eis given at the inlet and the static pressure is specified at the exit. 
Assumptions: steady and frictionless and the density of the fluid is constant.

Steps To Solve: Discritize, initialize, then calculate the new u using the A and b matrix maker and solver, then do the same for the pressure. 

=#
using Plots, FLOWMath 
include("SIMPLE1Funcs.jl")

# Givens
ρ = 1 #kg/s 
L = 2.0 #m
n = 4
dx = L/n #m
A_A = .5 #m^2
A_E = .1 #m^2
#Boundary Conditions
P0 = 10 #Pa 
Pe = 0 #Pa 
#Initial Conditions
mdot = 1.0 #kg/s 
#Under Relaxation Factors
αu = .8 #under relaxation factor for velocity
αp = .8 #under relaxation factor for pressure

#Discritize 
Pg = range(0,L,step=dx)
Vg = range(dx/2,L-dx/2,step=dx)

#Initial Conditions/guess
u = @. mdot/(ρ*Area(Vg)) #m/s, velocity at each velocity point
P = range(P0, Pe, length = n+1) #Pa, pressure at each pressure point

#Iterate until max(Apstar .* ustar - bstar) = 0
ustored, Pstored = iterator(u,P)

#prep for convergence study 
nalpha = [4,8,16,32,200,256]
prepalphas = [.8, .5, .2, .05, .01,.01]
alphafunc(xf) = linear(nalpha, prepalphas,xf)

#convergence study
function setup(n)
    global dx, Pg, Vg, u, P, αu, αp, mdot
    α = alphafunc(n)
    αp = α
    αu = α
    dx = L/n
    #Discritize 
    Pg = range(0,L,step=dx)
    Vg = range(dx/2,L-dx/2,step=dx)

    #Initial Conditions/guess
    u = @. mdot/(ρ*Area(Vg)) #m/s, velocity at each velocity point
    P = range(P0, Pe, length = n+1) #Pa, pressure at each pressure point
end
mdot4 = ρ * ustored[end][1] * Area(Vg[1])
setup(8)
println("n = 8")
ustored, Pstored = iterator(u,P)
mdot8 = ρ * ustored[end][1] * Area(Vg[1])
setup(16)
println("n = 16")
ustored, Pstored = iterator(u,P)
mdot16 = ρ * ustored[end][1] * Area(Vg[1])
setup(32)
αu = .05
αp = .05
println("n = 32")
ustored, Pstored = iterator(u,P)
mdot32 = ρ * ustored[end][1] * Area(Vg[1])
setup(200)
println("n = 200")
ustored, Pstored = iterator(u,P)
mdot200 = ρ * ustored[end][1] * Area(Vg[1])

order, Predicted = RichardsonExtrapolation([mdot4, mdot8, mdot16, mdot32],[4,8,16,32])

ConvgPlot = plot([4,8,16,32,200],[mdot4,mdot8,mdot16,mdot32,mdot200],label="mdot vs n",xlabel="n",ylabel="mdot (kg/s)",title="Convergence Study")
hline!([Predicted],label="Richardson Extrapolation Predicted mdot",style=:dash,color=:red)

#Analytical Pressure Solution Plot
#use bernoulli equation and ρ to solve for pressure at each point
mdotana = .44721 #kg/s
Pana = @. P0 - .5*ρ*mdotana^2 / ((ρ * Area(Pg))^2)

PPlot = plot(Pg,Pstored[end],label="Numerical Pressure",xlabel="Position (m)",ylabel="Pressure (Pa)",title="Pressure vs Position")
plot!(Pg,Pana,label="Analytical Pressure",style=:dash)



#=
ns = [4, 8, 16, 32, 64, 128, 256]
mdots = zero(Float64.(ns))
for (i, n) in enumerate(ns)
    global dx, Pg, Vg, u, P, αu, αp, mdot
    α = alphafunc(n)
    αp = α
    αu = α
    dx = L/n
    #Discritize 
    Pg = range(0,L,step=dx)
    Vg = range(dx/2,L-dx/2,step=dx)

    #Initial Conditions/guess
    u = @. mdot/(ρ*Area(Vg)) #m/s, velocity at each velocity point
    P = range(P0, Pe, length = n+1) #Pa, pressure at each pressure point
    un, pn = iterator(u,P)

    mdots[i] = ρ * un[end][1] * Area(Vg[1])
end
=#