#=
HW2FiniteDifFunctions.jl
February 1st, 2024
Jacob Child
Pseudocode: Define functions to use in HW2.jl to solve the 1D heat equation using the finite difference method. Currently planning to make it vector and matrix based, maybe eventually I can do structs.
=#

"""
Practice B: This function takes in bounds and the number of control volumes. It returns the node locations (length CV + 2), control surface locations (length CV + 1), and dx.
"""
function PracticeB(boundsf::Vector{Float64}, nf::Int)
    dx = (boundsf[2] - boundsf[1])/nf #Control volume width from splitting the bounds into n parts
    #Practice B puts nodes at the bounds and in the middle of the control volumes
    Interiors = boundsf[1] + dx/2:dx:boundsf[2] - dx/2 #interior node locations
    Nodes = [boundsf[1]; Interiors; boundsf[2]] #Node locations with the bounds
    CS = boundsf[1]:dx:boundsf[2] #Control surface locations (evenly spaced including the bounds)

    #calculate fw and fe for each node for the harmonic mean
    fw = zeros(length(Nodes))
    fe = zeros(length(Nodes))
    fw[1] = 0
    fe[1] = 1
    fw[end] = 1
    fe[end] = 0
    for i in 2:length(Nodes)-1
        fw[i] = (CS[i-1] - Nodes[i-1]) / (Nodes[i] - Nodes[i-1])
        fe[i] = (Nodes[i+1] - CS[i]) / (Nodes[i+1] - Nodes[i])
    end
    return Nodes, CS, dx #fw, fe, dx
end


"""
HMKw: Harmonic Mean Kw: This function takes in the thermal conductivity of the material in the control volume west of the node and the thermal conductivity of the primary CV, It also takes keyword arguments dx, and for dxw- and dxw+ that default to dx/2. It returns the harmonic mean of the thermal conductivities.
"""
function HMKw(kwf, kpf, dxf)
    dxw = dxf
    dxwminus = dxf/2
    dxwplus = dxf/2
    Kw = dxw*kwf*kpf / (kpf*dxwminus + kwf*dxwplus)
    #println("Kw is ", Kw)
    return Kw
end
function HMKsNew(kf,fwf, fef)
    kw = zeros(length(kf))
    ke = zeros(length(kf))
    for i in 1:length(kf)
        if i == 1
            ke[i] = 1 / ((1 - fef[i])/kf[i] + fef[i]/kf[i+1])
            kw[i] = 0 # this shouldn't be used 1 / ((1 - fwf[i])/kf[i] + fwf[i]/kf[i-1])
        elseif i == length(kf)
            ke[i] = 0 #This shouldn't be used 1 / ((1 - fef[i])/kf[i] + fef[i]/kf[i+1])
            kw[i] = 1 / ((1 - fwf[i])/kf[i] + fwf[i]/kf[i-1])
        else
            kw[i] = 1 / ((1 - fwf[i])/kf[i] + fwf[i]/kf[i-1])
            ke[i] = 1 / ((1 - fef[i])/kf[i] + fef[i]/kf[i+1])
        end
    end

    return kw, ke
end

function HMKw2ndNode(kwf, kpf, dxf)
    dxw = dxf
    dxwminus = 0
    dxwplus = dxf/2
    return dxw*kwf*kpf / (kpf*dxwminus + kwf*dxwplus)
end

function HMKwLastNode(kwf, kpf, dxf)
    dxw = dxf
    dxwminus = dxf/2
    dxwplus = 0
    return dxw*kwf*kpf / (kpf*dxwminus + kwf*dxwplus)
end

"""
HMKe: Harmonic Mean Ke: This function takes in the thermal conductivity of the material in the control volume east of the node and the thermal conductivity of the primary CV, It also takes keyword arguments dx, and for dxe- and dxe+ that default to dx/2. It returns the harmonic mean of the thermal conductivities.
"""
function HMKe(kef, kpf, dxf)
    dxe = dxf
    dxeminus = dxf/2
    dxeplus = dxf/2
    return dxe*kef*kpf / (kef*dxeminus + kpf*dxeplus)
end

function HMKeFirstNode(kef, kpf, dxf)
    dxe = dxf
    dxeminus = dxf/2
    dxeplus = 0
    return dxe*kef*kpf / (kef*dxeminus + kpf*dxeplus)
end

function HMKe2ndToLastNode(kef, kpf, dxf)
    dxe = dxf
    dxeminus = dxf/2
    dxeplus = 0
    return dxe*kef*kpf / (kef*dxeminus + kpf*dxeplus)
end

"""
MaterialMaker: This function takes in the vector of nodes, the thermal conductivities of the materials, and the location of the interface between the materials. It returns a vector of the thermal conductivities in the control volumes and a vector of the heat generation rates in the control volumes.
"""
function MaterialMaker(Nodesf, kf, qdotf, MatInterfaceLocf)
    k = zero(Nodesf) #Vector of thermal conductivities in the control volumes
    qdot = zero(Nodesf) #Vector of heat generation rates in the control volumes
    for (i, node) in enumerate(Nodesf)
        if node < MatInterfaceLocf
            k[i] = kf[1]
            qdot[i] = qdotf[1]
        else
            k[i] = kf[2]
            qdot[i] = qdotf[2]
        end
    end
    return k, qdot
end



function MaterialMaker5(Nodesf, kf, qdotf, MatInterfaceLocf)
    k = zero(Nodesf) #Vector of thermal conductivities in the control volumes
    qdot = zero(Nodesf) #Vector of heat generation rates in the control volumes
    for (i, node) in enumerate(Nodesf)
        if node < MatInterfaceLocf
            k[i] = kf[1]
            qdot[i] = qdotf[1]
        else
            k[i] = kf[2](node)
            qdot[i] = qdotf[2]
        end
    end
    return k, qdot
end





"""
FiniteDif: Finite Difference Function for 1D heat conduction with an adiabatic surface at the left boundary and a convection surface at the right boundary. The goal is to find the aw, ap, ae and b coefficients for each node and then solve the matrix to find the temperatures at each node. 
Inputs:
    kf: Vector of thermal conductivities in the materials (same length as Nodesf)
    qdotf: Vector of heat generation rates in the materials (same length as Nodesf)
    hf: Convection coefficient at the right boundary
    Tinf: Ambient temperature at the right boundary
    dxf: Width of the control volumes
    Nodesf: Vector of node locations
    CSf: Vector of control surface locations
    MatInterfaceLocf: Vector of the locations of the interfaces between materials
"""
function FiniteDif(kf, qdotf, Nodesf, hf, Tinf, dxf)
    
    #Initialize the vectors for the coefficients
    ap = zeros(length(Nodesf))
    ae = zeros(length(Nodesf))
    aw = zeros(length(Nodesf))
    b = zeros(length(Nodesf))
    #big for loop to go through the nodes
    #always compute the harmonic mean for ke and kw, and use MatInterfaceLoc to know when to switch materials and thus qdot and which k value is east and west etc
    for (i, _) in enumerate(Nodesf)
        #for the HMKe and HMKw function keyword arguments
        dx = dxf
        if i == 1 #Left boundary
            ap[i] = HMKeFirstNode(kf[i+1], kf[i], dxf) / dxf #!Why not 2*Ke?
            ae[i] = HMKeFirstNode(kf[i+1], kf[i], dxf) / dxf
            aw[i] = 0
            b[i] = 0
        elseif i == 2 #Node 2
            ap[i] = HMKe(kf[i+1], kf[i], dxf) / dxf + HMKw2ndNode(kf[1], kf[1], dxf) / dxf
            ae[i] = HMKe(kf[i+1], kf[i], dxf) / dxf
            aw[i] = HMKw2ndNode(kf[i-1], kf[i], dxf) / dxf
            b[i] = qdotf[1]*dxf
        elseif i == length(Nodesf) - 1 #Node 6
            ap[i] = HMKe2ndToLastNode(kf[i+1], kf[i], dxf) /dxf + HMKw(kf[i-1], kf[i], dxf) / dxf #!Why not 2*Ke?
            ae[i] = HMKe2ndToLastNode(kf[i+1], kf[i], dxf) / dxf
            aw[i] = HMKw(kf[i-1], kf[i], dxf) / dxf
            b[i] = 0
            #println("node6")
        elseif i == length(Nodesf) #Right boundary
            ap[i] = HMKwLastNode(kf[i-1], kf[i], dxf) / dxf + hf #!Why not 2*Kw?
            ae[i] = 0
            aw[i] = HMKwLastNode(kf[i-1], kf[i], dxf) / dxf
            b[i] = hf*Tinf
        else #Nodes 3-5, because I used MaterialMaker, I don't need to worry about where the conditions change
            ap[i] = HMKe(kf[i+1], kf[i], dxf) / dxf + HMKw(kf[i-1], kf[i], dxf) / dxf
            #println("At node $i, Ke is $(kf[i+1]) and Kw is $(kf[i-1]), and Kp is $(kf[i]).")
            ae[i] = HMKe(kf[i+1], kf[i], dxf) / dxf
            aw[i] = HMKw(kf[i-1], kf[i], dxf) / dxf
            b[i] = qdotf[i]*dxf
        end
    end
    #Make the A and b Matricies, solve for T 
    A = zeros(length(Nodesf), length(Nodesf))
    for i in 1:length(Nodesf)
        A[i, i] = ap[i]
        if i != length(Nodesf)
            A[i, i+1] = -ae[i]
        end
        if i != 1
            A[i, i-1] = -aw[i]
        end
    end
    T = A\b
    return T, A
end

"""
EasyRun: This function takes in the number of control volumes, the bounds, the thermal conductivities of the materials, the heat generation rates in the materials, the convection coefficient at the right boundary, and the ambient temperature at the right boundary. It returns the temperatures at each node.
"""
function EasyRun(Nf, Bf, kf, qdotf, hf, Tinf, MatInterfaceLocf)
    Nodes, CS, dx = PracticeB(Bf, Nf)
    ks, qdots = MaterialMaker(Nodes, kf, qdotf, MatInterfaceLocf)
    T, A = FiniteDif(ks, qdots, Nodes, hf, Tinf, dx)
    return T
end

function EasyRun5(Nf, Bf, kf, qdotf, hf, Tinf, MatInterfaceLocf)
    Nodes, CS, dx = PracticeB(Bf, Nf)
    ks, qdots = MaterialMaker5(Nodes, kf, qdotf, MatInterfaceLocf)
    T, A = FiniteDif(ks, qdots, Nodes, hf, Tinf, dx)
    return T
end

"""
ConvergenceStudy: This function takes in everything EasyRun5 takes in, plus a stopping criterion for the max temperature difference. It returns a vector of the max temperatures, a vector of the number of control volumes.
"""
function ConvergenceStudy(Nf, Bf, kf, qdotf, hf, Tinf, MatInterfaceLocf, StoppingCriterion)
    Tmax = 0
    Tmax2 = 1
    N = Nf
    Tmaxs = []
    Ns = []
    while abs(Tmax - Tmax2) > StoppingCriterion
        Tmax = Tmax2
        T = EasyRun5(N, Bf, kf, qdotf, hf, Tinf, MatInterfaceLocf)
        Tmax2 = maximum(T)
        push!(Tmaxs, Tmax2)
        push!(Ns, N)
        N += 1
    end
    return Tmaxs, Ns
end