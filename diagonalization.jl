module diagonalization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalize
#import potential
#import norm

function diagonalize(n,om,r,lambda,delta)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j for i in 0:n for j in -1:0]
 #println(vdiag)
 HMatrix=Array(Diagonal(vdiag))
 # non diagonal matrix elements
 #HMatrix[2,1]=1
 #HMatrix[1,2]=1
 #println(HMatrix[2,1])
 for i in 1:n
     HMatrix[2i+1,2i]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
     HMatrix[2i,2i+1]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
     HMatrix[2+2i,2i-1]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
     HMatrix[2i-1,2+2i]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
 end
 evals=eigvals(HMatrix)
 evecs=eigvecs(HMatrix)
 return [evals,HMatrix]
end

#test=diagonalize(50,1,0.1,0,1)
#println(test[1][1])
#println(test[1][2])
#println(test[2][6,3])
#println(test[2][3,6])

end