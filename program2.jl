using ITensors

function main()
    let
        
        # s = siteind("S=1/2")
        # Sx = op("Sx",s) 
        # @show Sz

        s = Index(2,"s")
        Sx = ITensor(s,prime(s))
        Sx[s[1],prime(s)[2]] = 0.5
        Sx[s[2],prime(s)[1]] = 0.5
        
        psi = ITensor(s)
        @show Sx

        psi[s[1]]=1
        psi[s[2]]=1
        psi = psi./norm(psi)
        println(psi)

        # D, U = eigen(Sx, s, prime(s))
        # println("D")
        # println(D)

        println("Sx psi")
        println(Sx*psi)

    end
end