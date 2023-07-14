using ITensors

# include("program3.jl")
# main()

function pr3_1()
    let
        
        # s = siteind("S=1/2")
        # Sx = op("Sx",s) 
        # @show Sz

        s = Index(2,"s")
        Sx = ITensor(s,prime(s))
        Sx[s[1],prime(s)[2]] = 0.5
        Sx[s[2],prime(s)[1]] = 0.5
        
        Sz = ITensor(s,prime(s))
        Sz[s[1],prime(s)[2]] = 0.5
        Sz[s[2],prime(s)[1]] = -0.5

        # psi = ITensor(s)
        # @show Sx

        # psi[s[1]]=1
        # psi[s[2]]=1
        # psi = psi./norm(psi)
        # println("psi")
        # println(psi)

        # # D, U = eigen(Sx, s, prime(s))
        # # println("D")
        # # println(D)

        # println("Sx psi")
        # println(Sx*psi)
        
        # cpsi = dag(prime(psi))

        # println("cpsi Sx psi")
        # println(cpsi*Sx*psi)


        s1 = Index(2,"s1")
        s2 = Index(2,"s2")

        psi = ITensor(s1,s2)
        psi[s1[1],s2[2]]=1.
        psi[s1[2],s2[1]]=-1.9 #for gs it was -1 but w can deviate for that and evolve to check if we reach gs for big beta

        psi = psi/sqrt(2.)
        psi = psi./norm(psi)
        println("psi")
        println(psi)

        Sz1 = ITensor(s1,prime(s1))
        Sz1[s1[1],prime(s1)[1]] = 0.5
        Sz1[s1[2],prime(s1)[2]] = -0.5
        Sz2 = ITensor(s2,prime(s2))
        Sz2[s2[1],prime(s2)[1]] = 0.5
        Sz2[s2[2],prime(s2)[2]] = -0.5

        Sp1 = ITensor(s1,prime(s1))
        Sp1[s1[1],prime(s1)[2]] = 1
        Sm1 = ITensor(s1,prime(s1))
        Sm1[s1[2],prime(s1)[1]] = 1

        println(Sp1)
        println(Sm1)

        Sp2 = ITensor(s2,prime(s2))
        Sp2[s2[1],prime(s2)[2]] = 1
        Sm2 = ITensor(s2,prime(s2))
        Sm2[s2[2],prime(s2)[1]] = 1

        println(Sp2)
        println(Sm2)

        H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2

        psiH = H*psi
        println("H psi")
        println(psiH)

        psibeta = exp(-100*H,ishermitian=true)*psi
        psibeta = psibeta./norm(psibeta)
        println("psibeta")
        println(psibeta)

        println("\n\n svd\n")
        U,S,V = svd(psibeta,s1,maxdim=1,cutoff=1.)   # compute SVD with (i,k) as row indices (indices of U)
        @show U
        @show S
        @show V
        println("psibeta")
        @show psibeta
        println("psibeta from svd")
        @show U*S*V

        println("\noverlap")
        newpsi = U*S*V
        @show dag(psibeta)*newpsi #check the overalp


    end
end