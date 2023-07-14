let
    N = 20
    sites = siteinds("S=1/2", N; conserve_qns=false)
    
    os = OpSum()
    for j in 1:N-1
        # os += 1.0,"Sz",j,"Sz",j+1
        os += 0.5,"S+",j,"S-",j+1
        os += 0.5,"S-",j,"S+",j+1
    end

    H = MPO(os, sites)

    psi0 = randomMPS(sites)
    sweeps = Sweeps(5)
    setmaxdim!(sweeps, 200)
    setcutoff!(sweeps, 1E-10)

    # E, psi = dmrg(H, psi0, sweeps, outputlevel = 0)
    # println("Ground state energy = $E")
    PH = ProjMPO(H)
    # for mpo it is position!
    

    for sweep in range(1,sweeps.nsweep)
        println(sweep)
        # for pos in range(1,N-1) # sic! WTF HOW DID IT MAKE IT WORK !!!!
        for pos in 1:N-1
            orthogonalize!(psi0,pos)  
            position!(PH,psi0,pos) #projmpo should be used here
            phi = psi0[pos]*psi0[pos+1]
            vals, vecs = eigsolve(PH,phi,1, ishermitian=true)
            energy = vals[1]
            phi::ITensor = vecs[1]
            U,S,V = svd(phi,inds(psi0[pos]), cutoff=sweeps.cutoff[sweep], maxdim=sweeps.maxdim[sweep])
            psi0[pos]=U
            psi0[pos+1] = S*V
            E = inner(psi0,apply(H,psi0))
            println("forward",pos,"       Energy: ",E)
        end
        # for pos in range(start=N-1,stop=1,step=-1)
        for pos in N-1:-1:1
            orthogonalize!(psi0,pos)  
            position!(PH,psi0,pos)
            phi = psi0[pos]*psi0[pos+1]
            vals, vecs = eigsolve(PH,phi,1)
            energy = vals[1]
            phi::ITensor = vecs[1]
            U,S,V = svd(phi,inds(psi0[pos]), cutoff=sweeps.cutoff[sweep], maxdim=sweeps.maxdim[sweep])
            psi0[pos+1]=V
            psi0[pos] = U*S
            E = inner(psi0,apply(H,psi0))
            println("backward",pos,"       Energy: ",E)
        end
        E = inner(psi0,apply(H,psi0))
        println("       Energy: ",E)
    end

end