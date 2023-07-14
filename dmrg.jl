using ITensors 
using KrylovKit: eigsolve

function sweeporder(n)
    result = []
    for i in 1:n-1
        push!(result, (i, :right))
    end
    for i in n-1:-1:1
        push!(result, (i, :left))
    end
    return result
end

let 
    N = 20
    sites = siteinds("S=1/2", N; conserve_qns=false)
    Δ = 1.0
    J = 1.0
  
    os = OpSum()
  
    for j in 1:(N - 1)
        os += Δ,"Sz",j,"Sz",j+1
        os += J/2,"S+",j,"S-",j+1
        os += J/2,"S-",j,"S+",j+1
    end

    H = MPO(os, sites)
    ψ = randomMPS(sites)
  
    sweeps = Sweeps(5)
    maxdim!(sweeps, 10, 20, 50, 100, 200)
    cutoff!(sweeps, 1E-10)

    PH = ProjMPO(H)
    for sw in 1:sweeps.nsweep
        for (pos, direction) in sweeporder(N)    
            
            orthogonalize!(ψ, pos)
            position!(PH, ψ, pos)

            ϕ = ψ[pos]*ψ[pos+1]

            vals, vecs = eigsolve(PH, ϕ, 1, ishermitian=true) # alternative for davidson

            energy = vals[1]
            ϕ::ITensor = vecs[1]

            U, S, V = svd(ϕ, inds(ψ[pos]), cutoff=sweeps.cutoff[sw], maxdim=sweeps.maxdim[sw])
            
            if direction == :right
                ψ[pos] = U
                ψ[pos+1] = S*V
            else
                ψ[pos] = U*S
                ψ[pos+1] = V
            end

            Epsi = apply(H, ψ)
            E = inner(ψ, Epsi)
            println("Sweep $sw, bond ($pos, $(pos+1)), direction $direction, energy = $E")
        end
    end
end