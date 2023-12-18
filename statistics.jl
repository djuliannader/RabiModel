module statistics
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
#import potential
#import norm

N=40
om=2
r=1
lambda=0.5
delta=0
nn=100
nu=2.0
chi=0.1

function analysisH(N,om,r,lambda,delta,nn,nu,chi)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu)
   eigvecsf=eigvecs(floquet)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta)
   evfvec = zeros(0)
   open("levels_output.dat","w") do io
   for i in 1:2*(N+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   evfvecs=sort(evfvec)
   for i in 1:length(evfvec)
     println(io,i," ",eigvs[1][i]," ",evfvecs[i])
   end
   end
   return "Done"
end

test= analysisH(N,om,r,lambda,delta,nn,nu,chi)
println("status: ",test)

end