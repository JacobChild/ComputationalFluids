#= 
TwoDSimpleFuncs.jl
Jacob Child
April 11th, 2024
Psuedocode: functions needed to assist in the SIMPLE CFD algorithm, also includes functions specific to the final 2D channel flow problem
=#

"""
uABMaker #TODO: finish the header
"""
function uABMaker(uf, vf, Pf; dxf = dx, dyf = dy, μf = μ, ρf = ρ, uinf = Vin, αuf = αu, wallflag = "internal")
    #uf is a 3 row matrix 
    #vf is a 2 row matrix
    #Pf is a Vector (1 row matrix), it should be the same row as the middle row in the uf matrix
    #!this code is made for a contant dxf (dx) and dyf (dy), but could be easily adjusted to be adaptable
    #East and West control surface values needed for the *internal velocity nodes* only
    De = μf/dxf
    Dw = μf/dxf
    Dn = μf/dyf
    Ds = μf/dyf
    
    #Initialize Ae, Aw, Ap, b, and d vectors to be put in the matrix later, use i,j I,J notation
    ai1J = zeros(size(uf,2)) #Ae
    aim1J = zeros(size(uf,2)) #Aw, aim1 = a(i-1)
    aiJ1 = zeros(size(uf,2)) #An
    aiJm1 = zeros(size(uf,2)) #As
    aiJ = zeros(size(uf,2)) #Ap
    BiJ = zeros(size(uf,2)) #b

    #Internal node for loop 
    J = 2 #ie the primary row we are looking at is the 2nd row of the 3 pulled in for u 
    internalend = length(aiJ)-1
    for i in 2:internalend #this leaves the boundaries out
        Fe = ρf / 2 * (uf[J,i+1] + uf[J, i])
        Fw = ρf / 2 * (uf[J,i] + uf[J, i-1])
        Fn = ρf / 2 * (vf[1,i+1] + vf[1,i]) #this row callout appears backwards from the equation. That is because I am starting in the top left, so j+1 is row 1 of v and j is row 2 of v
        Fs = ρf / 2 * (vf[2,i+1] + vf[2,i])
        
        ai1J[i] = (De*dyf + max(-Fe, 0)*dyf) #* αuf
        aim1J[i] = (Dw*dyf + max(Fw, 0)*dyf) #* αuf
        aiJ1[i] = (Dn*dxf + max(-Fn, 0)*dxf) #* αuf
        aiJm1[i] = (Ds*dxf + max(Fs, 0)*dxf) #* αuf
        aiJ[i] = (ai1J[i] + aim1J[i] + aiJ1[i] + aiJm1[i] +(Fe - Fw)*dyf + (Fn - Fs)*dxf) / αuf
        if wallflag == "lower"
            aiJm1[i] = 0
            aiJ[i] = (ai1J[i] + aim1J[i] + aiJ1[i] + aiJm1[i] +(Fe - Fw)*dyf + (Fn - Fs)*dxf + μf * dxf / (dyf/2)) / αuf

        elseif wallflag == "upper"
            aiJ1[i] = 0
            aiJ[i] = (ai1J[i] + aim1J[i] + aiJ1[i] + aiJm1[i] +(Fe - Fw)*dyf + (Fn - Fs)*dxf + μf * dxf / (dyf/2) ) / αuf
        end

        BiJ[i] = ((Pf[i-1] - Pf[i])*dyf + 0 + (1-αuf)/αuf * aiJ[i] * uf[J,i] * αuf)  #this is the u momentum equation with the underrelaxation factor added in
        
    end

    #println("Ap: ", Ap)
    #Boundary conditions/external nodes 
    #println(ai1J)
    #node 1
    
    ai1J[1] = 0
    aim1J[1] = 0
    aiJ1[1] = 0
    aiJm1[1] = 0
    aiJ[1] = 1
    BiJ[1] = uinf

    #node 5
    ai1J[end] = 0
    aim1J[end] = 1
    aiJ1[end] = 0
    aiJm1[end] = 0
    aiJ[end] = 1
    BiJ[end] = 0
    
    
    return ai1J, aim1J, aiJ1, aiJm1, aiJ, BiJ
end



function vABMaker(uf, vf, Pf; dxf = dx, dyf = dy, μf = μ, ρf = ρ, uinf = Vin, vinletf = 0, αvf = αv , wallflag = "internal")
    #uf is a 2 row matrix 
    #vf is a 3 row matrix amd row 2 is the row of interest I believe
    #Pf is a 2 row matrix, same location as uf
    #!this code is made for a contant dxf (dx) and dyf (dy), but could be easily adjusted to be adaptable
    #East and West control surface values needed for the *internal velocity nodes* only
    De = μf/dxf
    Dw = μf/dxf
    Dn = μf/dyf
    Ds = μf/dyf
    #Initialize Ae, Aw, Ap, b, and d vectors to be put in the matrix later, use i,j I,J notation
    ai1J = zeros(size(vf,2)) #Ae
    aim1J = zeros(size(vf,2)) #Aw, aim1 = a(i-1)
    aiJ1 = zeros(size(vf,2)) #An
    aiJm1 = zeros(size(vf,2)) #As
    aiJ = zeros(size(vf,2)) #Ap
    BiJ = zeros(size(vf,2)) #b

    #Internal node for loop 
    internalend = length(aiJ)-1 
    
    for i in 2:internalend #this leaves the boundaries out
        Fe = ρf / 2 * (uf[1,i] + uf[2, i])
        Fw = ρf / 2 * (uf[1,i-1] + uf[2, i-1])
        Fn = ρf / 2 * (vf[2,i] + vf[1,i]) #this row callout appears backwards, but should be correct
        Fs = ρf / 2 * (vf[3,i] + vf[2,i])
        
        if wallflag == "lower" || wallflag == "upper"
            ai1J[i] = 0
            aim1J[i] = 0
            aiJ1[i] = 0
            aiJm1[i] = 0
            aiJ[i] = 1
            BiJ[i] = 0
        else
        ai1J[i] = (De*dyf + max(-Fe, 0)*dyf) #* αvf
        aim1J[i] = (Dw*dyf + max(Fw, 0)*dyf) #* αvf
        aiJ1[i] = (Dn*dxf + max(-Fn, 0)*dxf) #* αvf
        aiJm1[i] = (Ds*dxf + max(Fs, 0)*dxf) #* αvf
        aiJ[i] = (ai1J[i] + aim1J[i] + aiJ1[i] + aiJm1[i] +(Fe - Fw)*dyf + (Fn - Fs)*dxf) / αvf
        
        BiJ[i] = (Pf[2,i-1] - Pf[1,i-1])*dxf  + 0 + (1-αvf)/αvf * aiJ[i] * vf[2,i] * αvf
        end
    end

    #println("Ap: ", Ap)
    #Boundary conditions/external nodes 
    #node 1
    ai1J[1] = 0
    aim1J[1] = 0
    aiJ1[1] = 0
    aiJm1[1] = 0
    aiJ[1] = 1
    BiJ[1] = vinletf #no y component of inlet velocity

    #node 6
    ai1J[end] = 0
    aim1J[end] = 1
    aiJ1[end] = 0
    aiJm1[end] = 0
    aiJ[end] = 1
    BiJ[end] = 0
    
    
    return ai1J, aim1J, aiJ1, aiJm1, aiJ, BiJ
end


function pPrimeABMaker(uapf, vapf, usnewf, vsnewf, Pf; dxf = dx, dyf = dy, αuf = αu, αvf = αv, αpf = αp, ρf = ρ, wallflag = "internal")
    #uapf is a single row vector, at the same row as P
    #vapf is a 2 row matrix, the row above the and below the P row of interest
    #usnewf is a single row vector at the same row as P
    #vsnewf is a 2 row matrix, the row above and below the P row of interest
    #!this code is made for a contant dxf (dx) and dyf (dy), but could be easily adjusted to be adaptable
    #we are solving at a single row of P

    #Initialize Ae, Aw, Ap, b, and d vectors to be put in the matrix later, use i,j I,J notation
    
    aIp1J = zeros(size(Pf,1)) #Ae
    aIm1J = zeros(size(Pf,1)) #Aw, aim1 = a(i-1)
    aIJp1 = zeros(size(Pf,1)) #An
    aIJm1 = zeros(size(Pf,1)) #As
    aIJ = zeros(size(Pf,1)) #Ap
    bpIJ = zeros(size(Pf,1)) #b

    #Internal node for loop
    internalend = size(aIJ,1)-1
    for i in 1:internalend #this leaves the boundaries out
        dip1J = dyf / uapf[i+1] #not multiplied by αuf as uapf already has it
        diJ = dyf / uapf[i] 
        dIjp1 = dxf / vapf[1,i+1] #not multiplied by αvf as vapf already has it
        dIj = dxf / vapf[2,i+1]
        
        if wallflag == "upper" && uapf[i] == 1 
            diJ = dyf / uapf[i] * αuf
            dIjp1 = dxf / vapf[1,i+1] * αvf
        end

        if wallflag == "upper" && uapf[i] != 1
            dIjp1 = dxf / vapf[1,i+1] * αvf
        end

        if wallflag == "internal" && uapf[i] == 1
            diJ = dyf / uapf[i] * αuf
        end

        if wallflag == "lower" && uapf[i] == 1
            diJ = dyf / uapf[i] * αuf
            dIj = dxf / vapf[2,i+1] * αvf
        end

        if wallflag == "lower" && uapf[i] != 1
            dIj = dxf / vapf[2,i+1] * αvf
        end

        aIp1J[i] = dip1J * ρf * dyf 
        aIm1J[i] = diJ * ρf * dyf 
        aIJp1[i] = dIjp1 * ρf * dxf 
        aIJm1[i] = dIj * ρf * dxf  #!off by a factor of 2, ie divide by 2?

        # #wall corrections
        # if wallflag == "lower"
        #     aIJm1[i] = 0 #these values are still included in aIJ, but are not used in the calculation overall for P
        # elseif wallflag == "upper"
        #     aIJp1[i] = 0
        # end

        aIJ[i] = aIp1J[i] + aIm1J[i] + aIJp1[i] + aIJm1[i]

        bpIJ[i] = -(ρf*usnewf[i+1]*dyf) + (ρf*usnewf[i]*dyf) - (ρf*vsnewf[1,i+1]*dxf) + (ρf*vsnewf[2,i+1]*dxf)
    
        
        #wall corrections
        if wallflag == "lower"
            aIJm1[i] = 0 #these values are still included in aIJ, but are not used in the calculation overall for P
        elseif wallflag == "upper"
            aIJp1[i] = 0
        end
    end

    #Boundary conditions/external nodes
    #node 5
    aIp1J[end] = 0
    aIm1J[end] = 0
    aIp1J[end-1] = 0
    aIJp1[end] = 0
    aIJm1[end] = 0
    aIJ[end] = 1
    bpIJ[end] = 0 #this is the exit condition I believe, pressure outlet

    #corner corrections
    if wallflag == "lower"
        aIm1J[1] = 0 #first node by the inlet has no aw
        aIp1J[end-1] = 0 #next to last node by the outlet has no ae
    elseif wallflag == "upper"
        aIm1J[1] = 0 #first node by the inlet has no aw
        aIp1J[end-1] = 0 #next to last node by the outlet has no ae
    end
    
    return aIp1J, aIm1J, aIJp1, aIJm1, aIJ, bpIJ

end


function iterator(u,v,P;)
    #Begin to solve, to start I will do 1 iteration at a time
    #preallocate us: east, west, north, south, B
    aues, auws, auns, auss, aups, Bus = [zero(u) for _ in 1:6]
    #preallocate vs: east, west, north, south, B
    aves, avws, avns, avss, avps, Bvs = [zero(v) for _ in 1:6]
    #Preallocate Ps: east, west, north, south, B 
    apes, apws, apns, apss, apps, Bps = [zero(P) for _ in 1:6]

    #upper wall 
    aues[1,:], auws[1,:], auns[1,:], auss[1,:], aups[1,:], Bus[1,:] = uABMaker(u[1:3,:], v[1:2,:], P[1,:], wallflag = "upper")

    #internal  for loop only 
    for i in axes(u[1:end-2,:],1)
        aues[i+1,:], auws[i+1,:], auns[i+1,:], auss[i+1,:], aups[i+1,:], Bus[i+1,:]  = uABMaker(u[i:i+2,:], v[i+1:i+2,:], P[i+1,:])
    end

    #bottom wall
    aues[end,:], auws[end,:], auns[end,:], auss[end,:], aups[end,:], Bus[end,:] = uABMaker(u[end-2:end,:], v[end-1:end,:], P[end,:], wallflag = "lower")

    Au = AMaker2D(aups, aues, auws, auns, auss, size(u))
    u1 = Au \ vec(Bus')
    u1 = reshape(u1, reverse(size(u)))' #all the transformations are required to get it looking correct
    #TODO: what is the proper way to do the mdot correction?
    mdotinapparent = sum(u1[:,1]) * W
    mdotoutappar = sum(u1[:,end]) * W
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

    Av = AMaker2D(avps, aves, avws, avns, avss, size(v))
    v1 = Av \ vec(Bvs');
    v1 = reshape(v1, reverse(size(v)))'

    #P calculations
    #upper wall
    apes[1,:], apws[1,:], apns[1,:], apss[1,:], apps[1,:], Bps[1,:] = pPrimeABMaker(aups[1,:],avps[1:2,:],u1[1,:],v1[1:2,:],P[1,:], wallflag = "upper")
    #internal  for loop only
    for i in axes(P[1:end-2,:],1)
    apes[i+1,:], apws[i+1,:], apns[i+1,:], apss[i+1,:], apps[i+1,:], Bps[i+1,:] = pPrimeABMaker(aups[i+1,:],avps[i+1:i+2,:],u1[i+1,:],v1[i+1:i+2,:],P[i+1,:])
    end
    #bottom wall 
    apes[end,:], apws[end,:], apns[end,:], apss[end,:], apps[end,:], Bps[end,:] = pPrimeABMaker(aups[end,:],avps[end-1:end,:],u1[end,:],v1[end-1:end,:],P[end,:], wallflag = "lower")
    println("here")
    Ap = AMaker2D(apps, apes, apws, apns, apss, size(P))
    P1 = Ap \ vec(Bps')
    P1 = reshape(P1, reverse(size(P)))'

    #Calculate the new P, u, v
    Pnew = P .+ αp .* P1
    diJ = dy ./ aups
    dIj = dx ./ avps
    unew = deepcopy(u1)
    vnew = deepcopy(v1)
    unew[:,2:end-1] = u1[:,2:end-1] .+ diJ[:,2:end-1] .* -diff(P1,dims=2)
    vnew[2:4,2:5] = v1[2:4,2:5] .+ dIj[2:4,2:5] .* diff(P1,dims=1)

    return unew, vnew, Pnew, Au, Av, Ap
end



function Converger(ustartf, vstartf, Pstartf;)
    #make all the keyword arguments global in this scope 

    #initialize variables to store 
    ustore, vstore, Pstore = [], [], []
    Aus, Avs, Aps = [], [], []
    Bus, Bvs, Bps = [], [], []
    
    ustore = push!(ustore, ustartf)
    vstore = push!(vstore, vstartf)
    Pstore = push!(Pstore, Pstartf) 
    u1, v1, P1, Au1, Av1, Ap1 = iterator(ustore[end], vstore[end], Pstore[end])
    push!(ustore, u1)
    push!(vstore, v1)
    push!(Pstore, P1)
    
    push!(Aus, Au1)
    push!(Avs, Av1)
    push!(Aps, Ap1)

    #iterate 
    #while !isapprox(1+maximum(Aus[end] * ustore[end-1] - Bus[end]), 1, atol = 1e-5) || !isapprox(1+maximum(Avs[end] * vstore[end-1] - Bvs[end]), 1, atol = 1e-5) || !isapprox(1+maximum(Aps[end] * Pstore[end-1] - Bps[end]), 1, atol = 1e-5)
    while maximum(ustore[end] - ustore[end-1]) > 1e-5 || maximum(vstore[end] - vstore[end-1]) > 1e-5 || maximum(Pstore[end] - Pstore[end-1]) > 1e-5

        u2, v2, P2, Au2, Av2, Ap2 = iterator(ustore[end], vstore[end], Pstore[end])
        push!(ustore, u2)
        push!(vstore, v2)
        push!(Pstore, P2)
        push!(Aus, Au2)
        push!(Avs, Av2)
        push!(Aps, Ap2)
        #to keep As and bs from growing huge
        popfirst!(Aus)
        popfirst!(Avs)
        popfirst!(Aps)
    end

    return ustore, vstore, Pstore
end

function AMaker2D(ap, ae, aw, an, as, sizer)
    ny, nx = sizer
    n = nx * ny
    A = zeros(n, n)
    for j in 1:ny
        for i in 1:nx
            k = (j-1)*nx + i
            A[k, k] = ap[j, i]  # Main diagonal
            if i != nx
                A[k, k+1] = -ae[j, i]  # East
            end
            if i != 1
                A[k, k-1] = -aw[j, i]  # West
            end
            if j != ny
                A[k, k+nx] = -an[j, i]  # North
            end
            if j != 1
                A[k, k-nx] = -as[j, i]  # South
            end
        end
    end
    return A
end



function RichardsonExtrapolation(ConMetricf, Nf)
    P = log((ConMetricf[3] - ConMetricf[2])/ (ConMetricf[4] - ConMetricf[3])) / log(2)
    println("The order of the simulation is $(round(P,digits=2)).")
    QFinal = ConMetricf[4] + (ConMetricf[3] - ConMetricf[4]) / (1 - 2^(round(P,digits=2)))
    println("The grid converged value is $(round(QFinal,digits=2)).")
    return P, QFinal
end