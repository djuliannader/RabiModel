module troterization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export troter
import dynamics
import diagonalization


function troter(Nmax::Int64,nn::Int64,r,om,lambda,delta,b,nu,eta,psi,flagt::Int64)
 pi=acos(-1)
 T=2*pi/nu
# b=nu*b
 dt=T/nn
#--- building the time independent Hamiltonian ----
 # diagonal matrix elements
 HMatrix = diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
 H1=HMatrix
 U0=exp(-im*H1*dt/2)
#--- builiding the time dependent Hamiltonian
 ddiag=[1 for i in 0:Nmax for j in -1:0]
 U=Diagonal(ddiag)
 if flagt==1
   utdiag=[b*(1/2)*(-1)^j for i in 0:Nmax for j in -1:0]
   Ut=Diagonal(utdiag)
 end
#
 if flagt==2
 vdiag=[0.0 for i in 0:Nmax for j in -1:0]
 Ut=Array(Diagonal(vdiag))
 for i in 1:Nmax
   Ut[2*i,2*i-1]=1.0*(b/2.0)*exp(im*psi)    # carrier
   Ut[2*i-1,2*i]=1.0*(b/2.0)*exp(-im*psi)   # carrier
 end
 end
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=(im/nu)*(-sin(i*nu*dt)+sin((i-1)*nu*dt))
  U1=exp(facu1*Ut)
  UM=U0*U1*U0
  U=UM*U
 end
 #println(H1)
 #println(U)
 return U
end

function troter2(Nmax::Int64,nn::Int64,r,om,gamma,omega,eta,psi)
 pi=acos(-1)
 T=2*pi/gamma
# b=nu*b
 dt=T/nn
#--- building the time dependent Hamiltonian RMP----
 ham=diagonalization.hamiltonian_rmp(Nmax,r,om,gamma,omega,eta,psi)
#--- building the time independent Hamiltonian ----
 h0=ham[1]
 U0=exp(-im*h0*dt/2)
#--- builiding the time dependent Hamiltonian
 Vt1=ham[2]
 Vt2=ham[3]
# ----
 ddiag=[1.0 + 0.0*im for i in 0:Nmax for j in -1:0]
 U=Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=1.0*(1.0/gamma)*(exp(-i*im*gamma*dt)-exp(-(i-1)*im*gamma*dt))
  facu2=1.0*(-1.0/gamma)*(exp(i*im*gamma*dt)-exp((i-1)*im*gamma*dt))
  #facu1=-im*dt
  #facu2=-im*dt
  opu=1.0*facu1*Vt1+1.0*facu2*Vt2
  U1=exp(1.0*opu)
  UM=U0*U1*U0
  U=UM*U
 end
 return U
end

function troter2_ct(Nmax::Int64,nn::Int64,r,om,gamma,omega,eta,psi)
 pi=acos(-1)
 T=2*pi/gamma
# b=nu*b
 dt=(T/nn)*im
#--- building the time dependent Hamiltonian RMP----
 ham=diagonalization.hamiltonian_rmp(Nmax,r,om,gamma,omega,eta,psi)
#--- building the time independent Hamiltonian ----
 h0=ham[1]
 U0=exp(-im*h0*dt/2)
#--- builiding the time dependent Hamiltonian
 Vt1=ham[2]
 Vt2=ham[3]
# ----
 ddiag=[1.0 + 0.0*im for i in 0:Nmax for j in -1:0]
 U=Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=1.0*(1.0/gamma)*(exp(-i*im*gamma*dt)-exp(-(i-1)*im*gamma*dt))
  facu2=1.0*(-1.0/gamma)*(exp(i*im*gamma*dt)-exp((i-1)*im*gamma*dt))
  #facu1=-im*dt
  #facu2=-im*dt
  opu=1.0*facu1*Vt1+1.0*facu2*Vt2
  U1=exp(opu)
  UM=U0*U1*U0
  U=UM*U
 end
 return U
end


function troter3(Nmax::Int64,nn::Int64,r,om,gamma,omega,eta,psi,no)
 pi=acos(-1)
 T=2*pi/gamma
# b=nu*b
 dt=(T/nn)
#--- building the time dependent Hamiltonian RMP----
 ham=diagonalization.hamiltonian_rmp2(Nmax,r,om,gamma,omega,eta,psi)
#--- building the time independent Hamiltonian ----
 h0=ham[1]
 U0=exp(-im*h0*dt/2)
 CCp=ham[4]
 CCm=ham[5]
#--- builiding the time dependent Hamiltonian
 a=ham[2]
 ad=ham[3]
# ----
 ddiag=[1.0 + 0.0*im for i in 0:Nmax for j in -1:0]
 U=Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
   opu=zeros(length(ddiag), length(ddiag))
   for j in 1:no
     facu1=(-im/(-im*gamma*j))*(exp(-i*im*gamma*dt*j)-exp(-(i-1)*im*gamma*dt*j))
     facu2=(-im/(im*gamma*j))*(exp(i*im*gamma*dt*j)-exp((i-1)*im*gamma*dt*j))
     #facu1=-im*dt
     #facu2=-im*dt
     opu=opu+(1/factorial(j))*((im*eta)^j*(CCp*(1.0*a*facu1+1.0*ad*facu2)^j)+(-im*eta)^j*(CCm*(1.0*ad*facu2+1.0*a*facu1)^j))
   end
  U1=exp(1.0*opu)  
  UM=U0*U1*U0
  U=UM*U
 end
 #k=1
 #hamtest = h0  + (1/factorial(k))*((im*eta)^k*(CCp*(1.0*a*1+1.0*ad*1)^k)+(-im*eta)^k*(CCm*(1.0*ad*1.0+1.0*a*1.0)^k))
 #evtest=eigvals(hamtest)
 #println("fraaag 3 ",evtest[1])
 return U
end




end
