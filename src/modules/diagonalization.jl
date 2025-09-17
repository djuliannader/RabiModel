module diagonalization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalize

function diagonalize(n,om,r,lambda,delta,eta,psi)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
 #println(vdiag)
 HMatrix=Array(Diagonal(vdiag))
 # non diagonal matrix elements
  for i in 1:n
     HMatrix[2*i+1,2*i]=1.0*(2/2)*eta*r^(1/2)*((delta+1)/2)*(i)^(1/2)*exp(-im*psi)   # Jaynes-Cummings
     HMatrix[2*i,2*i+1]=1.0*(2/2)*eta*r^(1/2)*((delta+1)/2)*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings
     HMatrix[2+2*i,2*i-1]=1.0*(2/2)*eta*r^(1/2)*((1-delta)/2)*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i-1,2+2*i]=1.0*(2/2)*eta*r^(1/2)*((1-delta)/2)*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i,2*i-1]=1.0*r*(lambda/2.0)*exp(im*psi)   # carrier
     HMatrix[2*i-1,2*i]=1.0*r*(lambda/2.0)*exp(-im*psi)   # carrier
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
     HMatrix[2*i+1,2*i]=1.0*(2/2)*eta*r^(1/2)*((delta+1)/2)*(i)^(1/2)*exp(-im*psi)   # Jaynes-Cummings
     HMatrix[2*i,2*i+1]=1.0*(2/2)*eta*r^(1/2)*((delta+1)/2)*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings
     HMatrix[2+2*i,2*i-1]=1.0*(2/2)*eta*r^(1/2)*((1-delta)/2)*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i-1,2+2*i]=1.0*(2/2)*eta*r^(1/2)*((1-delta)/2)*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings
     HMatrix[2*i,2*i-1]=1.0*r*(lambda/2.0)*exp(im*psi)   # carrier
     HMatrix[2*i-1,2*i]=1.0*r*(lambda/2.0)*exp(-im*psi)   # carrier
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
     Vt2[2*i+1,2*i]=-1.0*im*(omega/2)*eta*(i)^(1/2)*exp(-1.0*im*psi)   # Jaynes-Cummings         #3
     Vt1[2*i,2*i+1]=1.0*im*(omega/2)*eta*(i)^(1/2)*exp(im*psi)   # Jaynes-Cummings         #1
     Vt2[2+2*i,2*i-1]=1.0*im*(omega/2)*eta*(i)^(1/2)*exp(im*psi)   # Anti Jaynes-Cummings    #2
     Vt1[2*i-1,2+2*i]=-1.0*im*(omega/2)*eta*(i)^(1/2)*exp(-im*psi)   # Anti Jaynes-Cummings  #4
     H0[2*i,2*i-1]=1.0*(omega/2)*exp(im*psi)   # adds carrier
     H0[2*i-1,2*i]=1.0*(omega/2)*exp(-im*psi)   # adds carrier
 end
 return [H0,Vt1,Vt2]
end

function hamiltonian_rmp2(n::Int64,r,om,gamma,omega,eta,psi)
 # diagonal matrix elements
 vdiag=[(i)*om+r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
 vdiag2=[0.0+0.0*im for i in 0:n for j in -1:0]
 H0=Array(Diagonal(vdiag))
 a=Array(Diagonal(vdiag2))
 ad=Array(Diagonal(vdiag2))
 CCp=Array(Diagonal(vdiag2))
 CCm=Array(Diagonal(vdiag2))
 # non diagonal matrix elements
 for i in 1:(n)
     ad[2*i+1,2*i-1]=(i)^(1/2)*1.0   
     a[2*i-1,2*i+1]=(i)^(1/2)
     ad[2*i+2,2*i]=(i)^(1/2)
     a[2*i,2*i+2]=(i)^(1/2)   
     H0[2*i,2*i-1]=1.0*(omega/2)*exp(1.0*im*psi)   # add carrier
     H0[2*i-1,2*i]=1.0*(omega/2)*exp(-1.0*im*psi)   # add carrier
     CCp[2*i,2*i-1]=1.0*(omega/2)*exp(im*psi)   # only carrier
     CCm[2*i-1,2*i]=1.0*(omega/2)*exp(-1.0*im*psi)   # only carrier
 end
 return [H0,a,ad,CCp,CCm]
end


function sigmaz(n)
 # diagonal matrix elements
 vdiag=[(-1)^j+0*im for i in 0:n for j in -1:0]
 #println(vdiag)
 HMatrix=Array(Diagonal(vdiag))
 return HMatrix
end

function sigmax(n)
    vdiag=[0.0 for i in 0:n for j in -1:0]
    Sigx=Array(Diagonal(vdiag))
    for i in 1:n
      Sigx[2*i,2*i-1]=1.0
      Sigx[2*i-1,2*i]=1.0
    end
    return Sigx
end

end
