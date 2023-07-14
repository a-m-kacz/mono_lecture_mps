using ITensors
let

  


  i = Index(3)
  j = Index(5)
  k = Index(2)
  l = Index(7)

  A = ITensor(i,j,k)
  B = ITensor(j,l)

  # Set elements of A
  A[i=>1,j=>1,k=>1] = 11.1
  A[i=>2,j=>1,k=>2] = -21.2
  A[k=>1,i=>3,j=>1] = 31.1  # can provide Index values in any order
  # ...

  # Contract over shared index j
  C = A * B

  @show hasinds(C,i,k,l) # = true

  D = randomITensor(k,j,i) # ITensor with random elements

  # Add two ITensors
  # must have same set of indices
  # but can be in any order
  R = A + D

  nothing
  N = 10
  for j=1:N
    @show j
  end
end

# output

# hasinds(C, i, k, l) = true