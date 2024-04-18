#=
HW2.jl
February 1st, 2024
Jacob Child
Pseudocode: Solve the 1D heat equation using the finite difference method, first by hand, then using function calls etc with any number of nodes.
=#

using Plots

#Problem 3
#Part a with 7 nodes
ka = 15
kb = 60
qdota = 4*10^6
dx = .01
h = 1000 #W/(m^2*K)
T∞ = 300 #K
#Node 1
ap1 = 1
ae1 = 1
aw1 = 0
b1 = 0
#Node 2
ap2 = ka/dx + 2*ka/dx
ae2 = ka/dx
aw2 = 2*ka/dx
b2 = qdota*dx
#Node 3
ap3 = ka/dx + ka/dx
ae3 = ka/dx
aw3 = ka/dx
b3 = qdota*dx
#Node 4
kw4 = dx*ka*ka / (ka*dx/2 + ka*dx/2)
ke4 = dx*kb*ka / (kb*dx/2 + ka*dx/2)
ap4 = ke4/dx + kw4/dx
aw4 = kw4/dx 
ae4 = ke4/dx
b4 = qdota*dx
#Node 5
kw5 = dx*kb*ka / (kb*dx/2 + ka*dx/2)
ke5 = dx*kb*kb / (kb*dx/2 + kb*dx/2)
ap5 = ke5/dx + kw5/dx
aw5 = kw5/dx
ae5 = ke5/dx
b5 = 0
#Node 6
ap6 = 2*kb/dx + kb/dx
ae6 = 2*kb/dx
aw6 = kb/dx
b6 = 0
#Node 7 with convection boundary condition
ap7 = 2*kb/dx + h
ae7 = 0
aw7 = 2*kb/dx
b7 = h* T∞

#Part b make and solve the system of equations
A = [ap1 -ae1 0 0 0 0 0;
    -aw2 ap2 -ae2 0 0 0 0;
    0 -aw3 ap3 -ae3 0 0 0;
    0 0 -aw4 ap4 -ae4 0 0;
    0 0 0 -aw5 ap5 -ae5 0;
    0 0 0 0 -aw6 ap6 -ae6;
    0 0 0 0 0 -aw7 ap7]
b = [b1; b2; b3; b4; b5; b6; b7]
T = A\b

#Problem 4 (Same as Problem 3 a and b, but automated instead of by hand)
include("HW2FiniteDifFunctions.jl")
#Part a with 5 Control Volumes
bounds = [0,.05]
N = 5
#Material properties
ka = 15 #W/(m*K)
kb = 60 #W/(m*K)
qdota = 4*10^6 #W/m^3
h = 1000 #W/(m^2*K)
T∞ = 300 #K
MatInterfaceLoc = .03 #m
T4 = EasyRun(N, bounds, [ka,kb],[qdota,0], h, T∞, MatInterfaceLoc)
N2 = 20
T4b = EasyRun(N2, bounds, [ka,kb],[qdota,0], h, T∞, MatInterfaceLoc)
#Plot the difference between the two
x1,_,_ = PracticeB(bounds, N)
x2,CS2,dx2 = PracticeB(bounds, N2)
Plot4 = plot(x1, T4, label = "5 Control Volumes", xlabel = "X (m)", ylabel = "Temperature (K)", title = "Heat Distribution in 1D Wall", legend = :topleft)
plot!(x2, T4b, label = "20 Control Volumes", legend=:topright)

#Problem 5
#Same as before, but Kb is a function of x, so I will use MaterialMaker5 and EasyRun5 to solve it.
Kbfoo(xf) = 137*exp(25*xf-2) #Kb as a function of x
N2=20
T5 = EasyRun5(N, bounds, [ka,Kbfoo],[qdota,0], h, T∞, MatInterfaceLoc)
T5b = EasyRun5(N2, bounds, [ka,Kbfoo],[qdota,0], h, T∞, MatInterfaceLoc)
#Plot the difference between the two
Plot5 = plot(x1, T5, label = "5 Control Volumes", xlabel = "X (m)", ylabel = "Temperature (K)", title = "Heat Distribution in 1D Wall with Kb as a function of x", legend = :topleft)
plot!(x2, T5b, label = "20 Control Volumes", legend=:topright)

#Problem 5 part c, grid convergence study looking at the max temperature and iterating on N until Tmax is the same within .001
#Make a vector to store the Tmax values in so I can store the values for each N
#Make a while loop that runs while Tmax is not within .001 of the previous Tmax
#Tmaxs, Ns = ConvergenceStudy(N, bounds, [ka,Kbfoo],[qdota,0], h, T∞, MatInterfaceLoc, .001)
