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

function PecletNum(Nodesf, uf, rhof, gammaf, dxf)
    P = zero(Nodesf)
    F = zero(Nodesf)
    D = zero(Nodesf)
    for (i, node) in enumerate(Nodesf)
        F[i] = uf * rhof
        D[i] = gammaf/dxf[i]
        P[i] = F[i] / D[i]
    end
    return P, F, D
end

centfunc(pf) = 1-abs(pf)/2
upfunc(pf) = 1
hybfunc(pf) = max(0, 1-.5*abs(pf))
plfunc(pf) = max(0, (1-.1*abs(pf))^5)
SchemeFunctions = [centfunc, upfunc, hybfunc, plfunc]

function AbMakerHw5(Pf, Ff, Df, dxf, schemef)
    Ae = zero(Pf)
    Aw = zero(Pf)
    Ap = zero(Pf)
    b = zero(Pf)
    for (i, Pecf) in enumerate(Pf)
        A = SchemeFunctions[schemef](Pecf)
        if i == 1
            Ae[i] = 0
            Aw[i] = 0
            Ap[i] = 1
            b[i] = 1 #!initial condition given in this specific problem
        elseif i == length(Pf)
            Ae[i] = 0
            Aw[i] = 0
            Ap[i] = 1
            b[i] = 0 #!initial condition given in this specific problem
        else
            Ae[i] = Df[i] * SchemeFunctions[schemef](Pf[i]) + max(-Ff[i+1], 0)
            Aw[i] = Df[i] * SchemeFunctions[schemef](Pf[i]) + max(Ff[i-1], 0)
            Ap[i] = Ae[i] + Aw[i] #! no source term 
            b[i] = 0
        end
    end
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

