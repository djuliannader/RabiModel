module diagonalization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export troter
#import potential
#import norm

pi=3.1416
a=-pi
b=2.5*pi
c=-0.2*pi
om=2*pi
nn=1000

function troter(n::Int64,nn::Int64,a,b,c,om,pi)
 T=2*pi/om
 dt=T/nn
 Pauli3=[1 0; 0 -1]
 Pauli2=[0 -im; im 0]
 Pauli1=[0 1; 1 0]
 H1=a*Pauli3+c*Pauli1
 U0=exp(-im*H1*dt/2)
 U=Diagonal([1 for i in 1:n])
 for i in 1:nn
  facu1=(im*b/om)*(cos(i*om*dt)-cos((i-1)*om*dt))
  U1=exp(facu1*Pauli3)
  UM=U0*U1*U0
  U=UM*U
 end
 eigg=eigvals(U)
 return eigg
end

res=troter(2,nn,a,b,c,om,pi)
println(real(res[1]))


end