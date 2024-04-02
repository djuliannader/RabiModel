module diagonalization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalize
#import potential
#import norm

function diagonalize(n,om,r,lambda,delta,eta,psi)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
 #println(vdiag)
 HMatrix=Array(Diagonal(vdiag))
 # non diagonal matrix elements
  for i in 1:n
     HMatrix[2*i+1,2*i]=1.0*(2/2)*eta*((delta+1)/2)*(i)^(1/2)*exp(-im*psi)   # Jaynes-Cummings
     HMatrix[2*i,2*i+1]=1.0*(2/2)*eta*((delta+1)/2)*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings
     HMatrix[2+2*i,2*i-1]=1.0*(2/2)*eta*((1-delta)/2)*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i-1,2+2*i]=1.0*(2/2)*eta*((1-delta)/2)*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i,2*i-1]=1.0*(lambda/2)*exp(im*psi)   # carrier
     HMatrix[2*i-1,2*i]=1.0*(lambda/2)*exp(-im*psi)   # carrier
 end
 evals=eigvals(HMatrix)
 evecs=eigvecs(HMatrix)
 return [evals,HMatrix]
end

function hamiltonian(n,om,r,lambda,delta,eta,psi)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
 #println(vdiag)
 HMatrix=Array(Diagonal(vdiag))
 # non diagonal matrix elements
 for i in 1:n
     HMatrix[2*i+1,2*i]=1.0*(2/2)*eta*((delta+1)/2)*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings
     HMatrix[2*i,2*i+1]=1.0*(2/2)*eta*((delta+1)/2)*(i)^(1/2)*exp(-im*psi)   # Jaynes-Cummings
     HMatrix[2+2*i,2*i-1]=1.0*(2/2)*eta*((1-delta)/2)*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i-1,2+2*i]=1.0*(2/2)*eta*((1-delta)/2)*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i,2*i-1]=1.0*(lambda/2)*exp(im*psi)   # carrier
     HMatrix[2*i-1,2*i]=1.0*(lambda/2)*exp(-im*psi)   # carrier
 end
 return HMatrix
end



end