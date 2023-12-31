module troterization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export troter
import dynamics

#println("iniciando")
# --- parameters for troterization
#b=0.2 # correspons to chi in the article by Liu 2023
#om=5.0  # corresponds to nu in the article by Liu 2023
#nn=1000   # number of time intervals
# --- parameters or the time -independent hamiltonian
#Nmax=50
#om0=1.0  # fermionic frequency
#omc=2.0  # bosonic frequency
#lambda=0.5 # interaction
#delta=0.0 # JC - anti JC


function troter(Nmax::Int64,nn::Int64,om0,omc,lambda,delta,b,om)
 pi=acos(-1)
 T=2*pi/om
 b=om*b
 dt=T/nn
#--- building the time independent Hamiltonian ----
 # diagonal matrix elements
 vdiag=[(i)*omc+om0*(1/2)*(-1)^j for i in 0:Nmax for j in -1:0]
 HMatrix=Array(Diagonal(vdiag))
 for i in 1:Nmax
     HMatrix[2i+1,2i]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
     HMatrix[2i,2i+1]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
     HMatrix[2+2i,2i-1]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
     HMatrix[2i-1,2+2i]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
 end
 H1=HMatrix
 U0=exp(-im*H1*dt/2)
#--- builiding the time dependent Hamiltonian
 ddiag=[1 for i in 0:Nmax for j in -1:0]
 utdiag=[b*(1/2)*(-1)^j for i in 0:Nmax for j in -1:0]
 Ut=Diagonal(utdiag)
 U=Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=(im*b/om)*(-sin(i*om*dt)+sin((i-1)*om*dt))
  U1=exp(facu1*Ut)
  UM=U0*U1*U0
  U=UM*U
 end
 #println(H1)
 #println(U)
 return U
end

#Floquet=troter(Nmax,nn,om0,omc,lambda,delta,b,om)
#println("Floquet operator obtained")
#cs0=dynamics.initialcoherent(1.0,1.5,1.0,0.0,1.0,Nmax)
#println("initial coherent state built")
#mensaje1 = dynamics.survivalpt(cs0,Floquet,200.0,om)
#mensaje2 = dynamics.survivalp(cs0,200.0,1.0,Nmax,omc,om0,lambda,delta)
#println(mensaje2)


end
