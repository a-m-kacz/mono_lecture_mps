using ITensors
# using LinearAlgebra
using Plots
using KrylovKit: eigsolve 

function ex1()
    
    s = Index(2, "s")
    Sx = ITensor(s, prime(s))
    Sx[s[1], prime(s)[2]] = 0.5
    Sx[s[2], prime(s)[1]] = 0.5

    psi = ITensor(s)

    psi[s[1]] = 1.0
    psi[s[2]] = 1.0
    psi = psi./norm(psi)

    phi = Sx*psi
    println("phi: $phi")
    println("<psi|phi> = $(inner(psi,apply(Sx,psi)))")

end

function ex2()

    s1 = Index(2, "s1")
    s2 = Index(2, "s2")
    psi = ITensor(s1, s2)
    psi[s1(1), s2(2)] = 1.0/sqrt(2)
    psi[s1(2), s2(1)] = 1.0

    Sp1 = ITensor(s1, prime(s1))
    Sp1[s1(1), prime(s1)(2)] = 1.0
    Sp1[s1(2), prime(s1)(1)] = 0.0 

    Sm1 = ITensor(s1, prime(s1))
    Sm1[s1(1), prime(s1)(2)] = 0.0
    Sm1[s1(2), prime(s1)(1)] = 1.0

    Sp2 = ITensor(s2, prime(s2))
    Sp2[s2(1), prime(s2)(2)] = 1.0
    Sp2[s2(2), prime(s2)(1)] = 0.0

    Sm2 = ITensor(s2, prime(s2))
    Sm2[s2(1), prime(s2)(2)] = 0.0
    Sm2[s2(2), prime(s2)(1)] = 1.0

    Sz1 = ITensor(s1, prime(s1))
    Sz1[s1(1), prime(s1)(1)] = 0.5
    Sz1[s1(2), prime(s1)(2)] = -0.5

    Sz2 = ITensor(s2, prime(s2))
    Sz2[s2(1), prime(s2)(1)] = 0.5
    Sz2[s2(2), prime(s2)(2)] = -0.5

    H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2

    β = 100.
    expH = exp(-β*H, ishermitian=true)

    psibeta = expH*psi
    psibeta = noprime(psibeta)
    psibeta = psibeta./norm(psibeta)
    
    println("psibeta = ", psibeta)
    println("En = ", dag(prime(psibeta))*H*psibeta)

    ######### SVD #########
    U, S, V = svd(psibeta, s1, maxdim=1)
    newpsi = U*S*V
    overlap = dot(psibeta, newpsi)
    println("overlap for maxdim = 1: ", overlap)

    U, S, V = svd(psibeta, s1, maxdim=2)
    newpsi = U*S*V
    overlap = dot(psibeta, newpsi)
    println("overlap for maxdim = 2: ", overlap)
end

using LinearAlgebra
function ex3A()
    M = 
      [ 0.435839 0.223707 0.10;
        0.435839 0.223707 -0.10;
        0.223707 0.435839 0.10;
        0.223707 0.435839 -0.10 ]
    U, d, V = svd(M)
    
    dim = length(d)
    dim = 3
    Dtrunc = zeros((dim, dim))
    
    for i in 1:dim
        Dtrunc[i, i] = d[i]
    end

    Mtrunc = U*Dtrunc*(V')
    diff = norm((M - Mtrunc))
    println("|M - Mtrunc|^2 = ", diff^2)
    println("d = ", d)

end

function ex3B()
    s1 = Index(2, "s1")
    s2 = Index(2, "s2")

    sing = ITensor(s1, s2)
    prod = ITensor(s1, s2)

    sing[s1(1), s2(2)] = 1.0/sqrt(2)
    sing[s1(2), s2(1)] = -1.0/sqrt(2)

    prod[s1(1), s2(2)] = 1.0

    Ψ = ITensor(s1, s2)
    
    for mix in 0.0:0.1:1.0
        Ψ = (1-mix)*prod + mix*sing
        U, d, V = svd(Ψ, s1)
        if d[2, 2] ≈ 0.0
            entropy = -(d[1, 1]^2)*(log(ℯ, d[1, 1]^2))
        else
            entropy = -((d[1, 1]^2)*(log(ℯ, d[1, 1]^2)) + (d[2, 2]^2)*(log(ℯ, d[2, 2]^2)))
        end
        println("mix = ", mix, ",  vN entropy = ", entropy)
    end
end

function ex4()
    N = 50
    sites = siteinds("S=1/2", N, conserve_qns=false)

    h = 0.50
    
    # in Julia instead of ampo = AutoMPO(sites) we use OpSum() and after adding the operators 
    # we construct the hamiltonian as H = MPO(ampo, sites)
    
    os = OpSum()
    for j in 1:N-1
        os += -4.0,"Sz",j,"Sz",j+1
    end
    for j in 1:N
        os += -h*2.0,"Sx",j
    end

    H = MPO(os, sites)

    psi0 = randomMPS(sites)
    sweeps = Sweeps(5)
    setmaxdim!(sweeps, 20)
    setcutoff!(sweeps, 1E-10)

    E, psi = dmrg(H, psi0, sweeps, outputlevel = 0)
    println("Ground state energy = $E")

    bra = dag(prime(psi, "Site"))
    ket = psi
    Npos = floor(Int, N/2)
    Szjop = op(sites, "Sz", Npos)
    Sz_expect = bra[Npos]*Szjop*ket[Npos]
    println("Sz at site ", Npos, " = ", Sz_expect[1])
    
    # Now we want to calculate the expectation value of Sz at all sites
    
    # os = OpSum()
    # for j in 1:N
    #     os += "Sz",j
    # end
    # Szop = MPO(os, sites)
    
    Sz_avg = 0.
    for j in 1:N      
        Szjop = op(sites, "Sz", j)
        ket = psi
        bra = dag(prime(psi, "Site"))
        Sz_expect = bra[j]*Szjop*ket[j]
        println("Sz at site ", j, " = ", Sz_expect[1])
        Sz_avg +=Sz_expect[1]
    end
    println("average Sz = ", Sz_avg/N)

end

function ex4_1()
    N = 50
    sites = siteinds("S=1/2", N, conserve_qns=false)
    Sz_avg_all = Vector{Float64}()
    h_range = LinRange(0, 2, 50)
    for h in h_range        
        os = OpSum()
        for j in 1:N-1
            os += -4.0,"Sz",j,"Sz",j+1
        end
        for j in 1:N
            os += -h*2.0,"Sx",j
        end

        H = MPO(os, sites)

        psi0 = randomMPS(sites)
        sweeps = Sweeps(5)
        setmaxdim!(sweeps, 20)
        setcutoff!(sweeps, 1E-10)

        E, psi = dmrg(H, psi0, sweeps, outputlevel = 0)
        println("Ground state energy = $E")
    
        Sz_avg = 0.
        for j in 1:N      
            # needed to move the cente of orthogonality to j
            orthogonalize!(psi,j)  
            Szjop = op(sites, "Sz", j)
            ket = psi
            bra = dag(prime(psi, "Site"))
            Sz_expect = bra[j]*Szjop*ket[j]
            # println("Sz at site ", j, " = ", Sz_expect[1])
            Sz_avg +=abs(Sz_expect[1])
        end
        println("average Sz = ", Sz_avg/N)
        append!(Sz_avg_all, abs(Sz_avg/N))
    end
    plot(h_range, Sz_avg_all)
end

function ex5_1_trotter()
    N = 20
    cutoff = 1E-12
    tau = 0.1
    ttotal = 5.0

    # t = 2. dt = 0.1

    # Make an array of 'site' indices
    s = siteinds("S=1/2", N; conserve_qns=true)

    # Make gates (1,2),(2,3),(3,4),...
    ampo = OpSum()
    for j=1:N-1
        ampo += "Sz",j,"Sz",j+1
        ampo += 1/2,"S+",j,"S-",j+1
        ampo += 1/2,"S-",j,"S+",j+1
    end
    H = MPO(ampo,s)

    gates = ITensor[]
    for j in 1:(N - 1)
        s1 = s[j]
        s2 = s[j + 1]
        hj =
        op("Sz", s1) * op("Sz", s2) +
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2)
        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))

    # Initialize psi to be a product state (alternating up and down)
    psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

    # c = div(N, 2) # center site

    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    for t in 0.0:tau:ttotal
        # Sz = expect(psi, "Sz"; sites=c)
        # println("$t $Sz")
        E = inner(psi,Apply(H,psi))
        println("$t $E")

        t≈ttotal && break

        psi = apply(gates, psi; cutoff)
        normalize!(psi)
    end
    
end

function ex6_1()
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