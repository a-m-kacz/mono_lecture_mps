using ITensors

function tutorial3()
    let
    M = [0.435839 0.223707 0.10;
    0.435839 0.223707 -0.10;
    0.223707 0.435839 0.10;
    0.223707 0.435839 -0.10]

    print(M)

    dim = length(d)

    s1 = Index(2,"s1")
    s2 = Index(2,"s2")
    
    sing = ITensor(s1,s2)
    prod = ITensor(s1,s2)

    sing[s1[1],s2[2]] = 1/sqrt(2)
    sing[s1[2],s2[1]] = -1/sqrt(2)

    prod[s1[1],s2[2]]=1

    mixwf = ITensor(s1,s2)

    Sv = Array{Float64}(undef, 10)
    for n in 1:10
        println(n*0.1)
        mixwf = (1-n*0.1)*prod + n*0.1*sing
        U,S,V = svd(mixwf,s2)
        print(S)
        print(S.^2)
        S = S.^2
        Sv[n]=-1*(S[1,1]*log(S[1,1])+S[2,2]*log(S[2,2]))

    end




    return
    end
end