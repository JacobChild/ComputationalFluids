#=
HW5Funcs.jl
March 6th, 2024
Jacob Child
=#
"""
Practice A: This function takes in bounds and the number of nodes. It returns the node locations (length N), control surface locations (length N-1), and dx. Practice A does evenly spaced nodes, with the control surfaces in the middle of the nodes.
"""
function PracticeA(boundsf::Vector{Float64}, nf::Int)
    dx = (boundsf[2] - boundsf[1])/(nf-1) #Control volume width from splitting the bounds into n parts
    #Practice A puts nodes at the bounds and in the middle of the control volumes
    Nodes = range(boundsf[1], boundsf[2], length=nf) #Node locations with the bounds
    CS = range(boundsf[1]+dx/2, boundsf[2]-dx/2, length=nf-1) #Control surface locations (evenly spaced including the bounds)
    longdx = dx * ones(nf)
    longdx[1] = dx/2
    longdx[end] = dx/2
    return Nodes, CS, longdx
end

function HW6anbMaker(Nodesf, dxf, Tpstarf, Toldf; ρf= ρ, cf = c, kf = k, T∞f =  T∞, Tsurrf=Tsurr, ϵf = ϵ, hf = h, σf=σ, Tb=Tb, Δtf = Δt, bthickf = bthick, Tbf = Tb)
    Ae = zero(Nodesf)
    Aw = zero(Nodesf)
    Ap = zero(Nodesf)
    b = zero(Nodesf)
    #println(" $hf, $Tpstarf, $Toldf")
    Apof = @. ρf*cf*dxf/Δtf
    Sc = @. 2/bthickf * hf * T∞f + 2/bthickf *ϵf*σf*Tsurrf^4 + 6*ϵf*σf/bthickf * Tpstarf^4
    Sp = @. -2/bthickf * hf - 8*ϵf*σf/bthickf * Tpstarf^3
    for (i, node) in enumerate(Nodesf)
        if i == 1
            Ae[i] = 0
            Aw[i] = 0
            Ap[i] = 1
            b[i] = Tbf #!initial condition given in this specific problem
        elseif i == length(Nodesf)
            Ae[i] = 0
            Aw[i] = kf/(2*dxf[i]) #my little dxf is still standard size
            Ap[i] = Apof[i] + Aw[i] - Sp[i]*dxf[i]
            b[i] = Sc[i]*dxf[i] + Apof[i]*Toldf[i] #!source specific to this problem
        else
            Ae[i] = kf/dxf[i-1]
            Aw[i] = kf/dxf[i+1]
            Ap[i] = Ae[i] + Aw[i] + Apof[i] - Sp[i]*dxf[i] #! no source term 
            b[i] = Sc[i]*dxf[i] + Apof[i]*Toldf[i] #!source specific to this problem
        end
    end
    #println("Tb is $Tbf")
    return Ae, Aw, Ap, b
end


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

function TimeStepper(Nodesf, dxf, Tinitf, Δtf, tmaxf, T∞f, Tsurrf, bthickf, hf, ϵf, σf, Tb, ρf, cf, kf)
    #rename givens 
    ρ = ρf
    c = cf
    k = kf
    T∞ = T∞f
    Tsurr = Tsurrf
    ϵ = ϵf
    σ = σf
    bthick = bthickf
    h = hf
    # Initialization
    T = ones(length(Nodesf)) * Tinitf
    Tstar = ones(length(Nodesf)) * Tinitf
    Told = ones(length(Nodesf)) * Tinitf
    # Create an array to store the temperature vectors
    T_all = []
    # Time loop
    t = 0
    while t < tmaxf
        # Iterative convergence solver for each time step
        while true
            #println("iterating at time $t")
            Tstar = copy(T)
            Ae, Aw, Ap, b = HW6anbMaker(Nodesf, dxf, Tstar, Told; ρf= ρ, cf = c, kf = k, T∞f =  T∞, Tsurrf=Tsurr, ϵf = ϵ, hf = h, σf=σ, Tb=Tb, Δtf = Δtf, bthickf = bthickf, Tbf = Tb)
            T = TDMASolver(Ap, Ae, Aw, b)
            if maximum(abs.(Tstar - T)) <= 10e-4
                break
            end
        end
        # Save the temperature vector at this time step
        push!(T_all, copy(T))
        Told = copy(T)
        t += Δtf
    end
    return T_all
end