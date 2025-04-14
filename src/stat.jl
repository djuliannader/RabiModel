module stat
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
import wigner_eig
using QuantumOptics

function analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   #println("flag",eigvs[1][1],eigvs[1][2])
   println("See file DensityOfStates_output.dat for Density of states (DOS) ")
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
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
   #open("DensityOfStates_output2.dat","w") do ioa
   # for i in 1:trunc(Int64,length(eigvs[1]))
   #   println(ioa,eigvs[1][i]/r)
   # end
   #end
   end
   end
   return "Done"
end


function parameter_r(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
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

 function orderinfvec(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
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

function orderinfvec3(Nmax,om,r,gamma,omega,eta,psi,nn,no)
   floquet=troterization.troter3(Nmax,nn,r,om,gamma,omega,eta,psi,no)
   hampart=diagonalization.hamiltonian_rmp2(Nmax,r,om,gamma,omega,eta,psi)
   ham=hampart[1]
   CCp=hampart[4]
   CCm=hampart[5]
   a = hampart[2]
   ad= hampart[3]
   #for j in 1:no
     j=1
     ham = ham  + (1/factorial(j))*((im*eta)^j*(CCp*(1.0*a*1+1.0*ad*1)^j)+(-im*eta)^j*(CCm*(1.0*ad*1.0+1.0*a*1.0)^j))
   #end
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
   #println(evford[1])   #  descomentar para ver como se ordenan las quasienergias
   #println("flaaag 2 lala", evfvec[evford[1]])   #  descomentar para ver como se ordenan las quasienergias
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

  function quasienergies(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
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

function purity_rmp(Nmax,om,r,omega,gamma,eta,psi,Nf)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   floquet=troterization.troter2(Nmax,Nf,r,om,gamma,omega,eta,psi)
   vecordered = orderinfvec2(Nmax,om,r,gamma,omega,eta,psi,Nf)
   evs=eigvecs(floquet)
   println("Purity of stationary states")
   purity=[]
   entropy=[]
   for k in 1:50
     listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
     phi = wigner_eig.buildingstate(listvec,Nmax)
     rho = dm(phi)
     rhoqp = ptrace(rho,2)
     pur = tr(rhoqp*rhoqp)
     ent = entropy_vn(rhoqp)
     append!(purity,pur)
     append!(entropy,ent)
    end
   return [purity,entropy] 
 end


function purity(Nmax,om,r,lambda,delta,eta,psi)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
   evs=eigvecs(ham)
   println("Purity of stationary states")
   purity=[]
   entropy=[]
   for k in 1:50
     listvec=[evs[i,k] for i in 1:2*Nmax]
     phi = wigner_eig.buildingstate(listvec,Nmax)
     rho = dm(phi)
     rhoqp = ptrace(rho,2)
     pur = tr(rhoqp*rhoqp)
     ent = entropy_vn(rhoqp)
     append!(purity,pur)
     append!(entropy,ent)
    end
    evv=eigvals(ham)
   return [purity,entropy,evv] 
 end

end