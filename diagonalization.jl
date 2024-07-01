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

function hamiltonian_rmp(n::Int64,r,om,gamma,omega,eta,psi)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
 vdiag2=[0.0+0.0*im for i in 0:n for j in -1:0]
 #println(vdiag)
 H0=Array(Diagonal(vdiag))
 Vt1=Array(Diagonal(vdiag2))
 Vt2=Array(Diagonal(vdiag2))
 # non diagonal matrix elements
 for i in 1:n
     Vt2[2*i+1,2*i]=im*(omega/2)*eta*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings
     Vt1[2*i,2*i+1]=-im*(omega/2)*eta*(i)^(1/2)*exp(-im*psi)   # Jaynes-Cummings
     Vt2[2+2*i,2*i-1]=im*(omega/2)*eta*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings
     Vt1[2*i-1,2+2*i]=-im*(omega/2)*eta*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings
     H0[2*i,2*i-1]=1.0*(omega/2)*exp(im*psi)   # carrier
     H0[2*i-1,2*i]=1.0*(omega/2)*exp(-im*psi)   # carrier
 end
 return [H0,Vt1,Vt2]
end


end