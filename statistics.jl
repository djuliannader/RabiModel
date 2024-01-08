module statistics
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
#import potential
#import norm


function analysisH(N,om,r,lambda,delta,nn,nu,chi)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta)
   open("DensityOfStates_output.dat","w") do ioa
   for i in 1:trunc(Int64,length(eigvs[1])/2-1)
      println(ioa,(eigvs[1][2*i+1]+eigvs[1][2*i])/2," ",1.0/(eigvs[1][2*i+1]-eigvs[1][2*i]))
   end
   end
   println("See file DensityOfStates_output.dat for Density of states (DOS) ")
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


function parameter_r(N,om,r,lambda,delta,nn,nu,chi)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu)
   eigvalsf=eigvals(floquet)
   quasienergies=zeros(0)
   for i in 1:length(eigvalsf)
      fase=atan(imag(eigvalsf[i])/real(eigvalsf[i]))
      append!(quasienergies, pi-fase )
   end
   quasienergies=sort(quasienergies)
   #println(quasienergies)
   sn=[quasienergies[n+1]-quasienergies[n] for n in 1:(length(quasienergies)-1)]
   #sn=[quasienergies[2*n+1]-quasienergies[2*n] for n in 1:trunc(Int64,length(quasienergies)/2-1)]
   rn=zeros(0)
   for n in 1:length(sn)-1
       if sn[n]>sn[n+1]
         append!(rn, sn[n+1]/sn[n])
       end
       if sn[n]<sn[n+1]
         append!(rn, sn[n]/sn[n+1])
       end
   end
   r=sum(rn)/length(rn)
   return r
   end




end