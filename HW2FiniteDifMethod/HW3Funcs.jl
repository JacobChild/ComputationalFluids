#=
HW3Funcs.jl
=#


"""
Practice B: This function takes in bounds and the number of control volumes. It returns the node locations (length CV + 2), control surface locations (length CV + 1), and dx.
"""
function PracticeB(boundsf::Vector{Float64}, nf::Int)
    dx = (boundsf[2] - boundsf[1])/nf #Control volume width from splitting the bounds into n parts
    #Practice B puts nodes at the bounds and in the middle of the control volumes
    Interiors = range(boundsf[1] + dx/2, boundsf[2] - dx/2, length=nf) #interior node locations
    #! The above line is different than the original, but it is more robust and less error prone
    #! original: Interiors = boundsf[1] + dx/2:dx:boundsf[2] - dx/2 #interior node locations
    Nodes = [boundsf[1]; Interiors; boundsf[2]] #Node locations with the bounds
    #CS = boundsf[1]:dx:boundsf[2] #Control surface locations (evenly spaced including the bounds)
    CS = range(boundsf[1], boundsf[2], length(Nodes)-1) #? The above causes issues at weird points like N = 44

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
    return Nodes, fw, fe, dx
end


function HMKNew(kf,fwf, fef)
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
function FiniteDif(kwf, kef, qdotf, Nodesf, hf, Tinf, dxf)
    
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
            ap[i] = 2*kef[i] / dxf
            ae[i] = 2*kef[i] / dxf
            aw[i] = 0
            b[i] = 0
        elseif i == 2 #Node 2
            ap[i] = kef[i] / dxf + 2*kwf[i] / dxf
            ae[i] = kef[i] / dxf
            aw[i] = 2*kwf[i] / dxf
            b[i] = qdotf[1]*dxf
        elseif i == length(Nodesf) - 1 #Node 6
            ap[i] = 2*kef[i] /dxf + kwf[i] / dxf 
            ae[i] = 2*kef[i] / dxf
            aw[i] = kwf[i] / dxf
            b[i] = 0
            
        elseif i == length(Nodesf) #Right boundary
            ap[i] = 2*kwf[i] / dxf + hf
            ae[i] = 0
            aw[i] = 2*kwf[i] / dxf
            b[i] = hf*Tinf
        else #Nodes 3-5, because I used MaterialMaker, I don't need to worry about where the conditions change
            ap[i] = kef[i] / dxf + kwf[i] / dxf
            #println("At node $i, Ke is $(kf[i+1]) and Kw is $(kf[i-1]), and Kp is $(kf[i]).")
            ae[i] = kef[i] / dxf
            aw[i] = kwf[i] / dxf
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
    return T, A, b
end


function EasyRun5(bounds, N, ka, kb, qdota, h, T∞, MatInterfaceLoc)
    x,fw, fe, dx = PracticeB(bounds, N)
    #safety check
    if length(x) -2 != N
        throw("PracticeB did not return the correct number of nodes. N = $N")
    end
    ks, qdot = MaterialMaker5(x,[ka,kb],[qdota, 0],MatInterfaceLoc)
    #safety check
    if length(ks) -2 != N
        throw("MaterialMaker5 did not return the correct number of materials. N = $N")
    end
    kw,ke = HMKNew(ks,fw,fe)
    NewKw = []
    NewKe = []
    for i in 1:length(ks)
        if i == 1
            push!(NewKw, "Kwcatch")
            push!(NewKe, ke[i])
        elseif i == length(ks)
            push!(NewKw, kw[i])
            push!(NewKe, "KeCatch")
        else
            push!(NewKw, kw[i])
            push!(NewKe, ke[i])
        end 
    end
    T, A, b = FiniteDif(NewKw, NewKe, qdot, x, h, T∞, dx)
    return T, A, b
end


###
###
#Functions for HW3 Problems 2-3 TDMA

function TDMASolver(a,b,c,d)
    n = length(d)
    P = zero(Float64.(d))
    Q = zero(Float64.(d))
    Φ = zero(Float64.(d))
    for i in 1:n
        if i == 1
            P[i] = b[i] / a[i]
            Q[i] = d[i] / a[i]
        else
            P[i] = b[i] / (a[i] - c[i]*P[i-1])
            Q[i] = (d[i] + c[i]*Q[i-1]) / (a[i] - c[i]*P[i-1])
        end
    end
    Φ[n] = Q[n]
    ns = 1:n
    for i in reverse(ns[1:end-1])
        Φ[i] = P[i]*Φ[i+1] + Q[i]
    end
    return Φ
end

"""
abmaker: function to make the ae, ap, aw, and b vectors given a Tpstar, K values, dx, and all the given constants
"""
function abmaker(Tpstar, Kf, dxf, Df, emf, sigmaf, hf, Tbasef, Tinff, Tsurrf)
    n = length(Kf)
    ap = zero(Float64.(Kf))
    ae = zero(Float64.(Kf))
    aw = zero(Float64.(Kf))
    b = zero(Float64.(Kf))
    #Calculations
    Sc = @. 4/Df *hf*Tinff + 4*emf*sigmaf/Df * Tsurrf^4 + 12 * emf * sigmaf * Tpstar^4 / Df
    Sp = @. -4 * hf / Df - 16 * emf * sigmaf * Tpstar^3 / Df
    for i in 1:n
        if i == 1
            aw[i] = 0
            ae[i] = 0
            ap[i] = ae[i] + aw[i] + 1
            b[i] = Tbasef
        elseif i == 2
            aw[i] = 2*Kf[i] / dxf
            ae[i] = Kf[i] / dxf
            ap[i] = ae[i] + aw[i] - Sp[i] * dxf 
            b[i] = Sc[i] * dxf
        elseif i == n-1
            aw[i] = Kf[i] / dxf
            ae[i] = 2*Kf[i] / dxf
            ap[i] = ae[i] + aw[i] - Sp[i] * dxf
            b[i] = Sc[i] * dxf
        elseif i == n
            #adiabatic tip
            aw[i] = 1
            ae[i] = 0
            ap[i] = ae[i] + aw[i]
            b[i] = 0 #? is this right?
        else
            aw[i] = Kf[i] / dxf
            ae[i] = Kf[i] / dxf
            ap[i] = ae[i] + aw[i] - Sp[i] * dxf
            b[i] = Sc[i] * dxf
        end
    end
    return aw, ap, ae, b
end
