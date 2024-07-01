module stat
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization

function analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   #println("flag",eigvs[1][1],eigvs[1][2])
   println("See file DensityOfStates_output.dat for Density of states (DOS) ")
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
   eigvecsf=eigvecs(floquet)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
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
   open("DensityOfStates_output.dat","w") do ioa
   if abs(lambda)>0.0
     for i in 1:trunc(Int64,length(eigvs[1])/2-1)
       println(ioa,(eigvs[1][2*i]+eigvs[1][2*i-1])/2," ",1.0/(eigvs[1][2*i]-eigvs[1][2*i-1])," ",(evfvecs[2*i]+evfvecs[2*i-1])/2," ",1.0/(evfvecs[2*i]-evfvecs[2*i-1]))
     end
   else
     for i in 1:trunc(Int64,length(eigvs[1])/2-1)
       println(ioa,(eigvs[1][2*i+1]+eigvs[1][2*i-1])/2," ",1.0/(eigvs[1][2*i+1]-eigvs[1][2*i-1])," ",(evfvecs[2*i+1]+evfvecs[2*i-1])/2," ",1.0/(evfvecs[2*i+1]-evfvecs[2*i-1]))
     end
   end
   end
   end
   return "Done"
end


function parameter_r(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
   eigvalsf=eigvals(floquet)
   qs=zeros(0)
   for i in 1:length(eigvalsf)
      fase=mod(real(log(eigvalsf[i])/(-im)),2*pi)
      append!(qs,fase)
   end
   qs=sort(qs)
   sn=[qs[n+1]-qs[n] for n in 1:(length(qs)-1)]
   fac=sum(sn)/length(sn)
   sn=(1/fac)*(sn)
   #println(sn)
   rn=zeros(0)
   open("quasienergies_spacing.dat","w") do io
     for i in 1:length(sn)
       println(io,sn[i])
     end
   end
   for n in 1:length(sn)-1
       if sn[n]>sn[n+1]
         append!(rn, sn[n+1]/sn[n])
       end
       if sn[n]<sn[n+1]
         append!(rn, sn[n]/sn[n+1])
       end
     #println(io," ",sn[n]) 
   end
   #end
   rpar=sum(rn)/length(rn)
   return rpar
   end

 function orderinfvec(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
   eigvecsf=eigvecs(floquet)
   evfvec = zeros(0)
   for i in 1:2*(N+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   evford=sortperm(evfvec)
   #println(evford)   #  descomentar para ver como se ordenan las quasienergias
   return evford
  end

function orderinfvec2(Nmax,om,r,gamma,omega,eta,psi,nn)
   floquet=troterization.troter2(Nmax,nn,r,om,gamma,omega,eta,psi)
   hampart=diagonalization.hamiltonian_rmp(Nmax,r,om,gamma,omega,eta,psi)
   ham=hampart[1]+hampart[2]+hampart[3]
   eigvecsf=eigvecs(floquet)
   evfvec = zeros(0)
   for i in 1:2*(Nmax+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(Nmax+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   evford=sortperm(evfvec)
   #println(evford)   #  descomentar para ver como se ordenan las quasienergias
   return evford
  end

function dos_rmp(Nmax,om,r,gamma,omega,eta,psi,nn)
   floquet=troterization.troter2(Nmax,nn,r,om,gamma,omega,eta,psi)
   hampart=diagonalization.hamiltonian_rmp(Nmax,r,om,gamma,omega,eta,psi)
   ham=hampart[1]+hampart[2]+hampart[3]
   eigvecsf=eigvecs(floquet)
   evfvec = zeros(0)
   for i in 1:2*(Nmax+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(Nmax+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   println(sort(evfvec))
   #evford=sortperm(evfvec)
   #println(evford)   #  descomentar para ver como se ordenan las quasienergias
   return "done"
  end

  function quasienergies(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
   #eigvalsf=eigvals(floquet)
   #qs=zeros(0)
   #for i in 1:length(eigvalsf)
   #   fase=mod(real(log(eigvalsf[i])/(-im)),2*pi)
   #   append!(qs,fase)
   #end
   #qs=sort(qs)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
   eigvecsf=eigvecs(floquet)
   evfvec = zeros(0)
   for i in 1:2*(N+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   ev_sorted=sort(evfvec)
   return ev_sorted
  end


 function orderinfvecqs(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
   eigvalsf=eigvals(floquet)
   qs=zeros(0)
   for i in 1:length(eigvalsf)
      fase=mod(real(log(eigvalsf[i])/(-im)),2*pi)
      append!(qs,fase)
   end
   qs=sort(qs)
   #ham=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
   #eigvecsf=eigvecs(floquet)
   #evfvec = zeros(0)
   #for i in 1:2*(N+1)
   #  eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
   #  eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
   #  for k in 1:length(eigvecf)
   #    eigvecf_t[1,k]=conj(eigvecf[k])
   #  end
   #  evf=eigvecf_t*ham*eigvecf
   #  append!(evfvec, real(evf) )
   #end
   evfsorted=sortperm(qs)
   return evfsorted
  end

end