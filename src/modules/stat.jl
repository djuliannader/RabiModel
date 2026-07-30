module stat
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("diagonalization.jl")
include("troterization.jl")
using .diagonalization
using .troterization



function analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   #println("flag",eigvs[1][1],eigvs[1][2])
   println("See file output/DensityOfStates_output.dat for Density of states (DOS) ")
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
   eigvecsf=eigvecs(floquet)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
   evfvec = zeros(0)
   open("output/levels_output.dat","w") do io
   for i in 1:2*(N+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf))
   end
   evfvecs=sort(evfvec)
   for i in 1:length(evfvec)
     println(io,i," ",eigvs[1][i]," ",evfvecs[i])
   end
   open("output/DensityOfStates_output.dat","w") do ioa
   if abs(lambda)>0.0
     for i in 1:trunc(Int64,length(eigvs[1])/2-1)
       println(ioa,(eigvs[1][2*i]+eigvs[1][2*i-1])/2," ",1.0/(eigvs[1][2*i]-eigvs[1][2*i-1])," ",(evfvecs[2*i]+evfvecs[2*i-1])/2," ",1.0/(evfvecs[2*i]-evfvecs[2*i-1]))
     end
   else
       #println("heeere")
       for i in 1:trunc(Int64,length(eigvs[1])/2-1)
        eav =  (eigvs[1][i+1]+eigvs[1][i])/(2*r)
       if eav > -0.5
           println(ioa,(eigvs[1][i+1]+eigvs[1][i])/2," ",1.0/(eigvs[1][i+1]-eigvs[1][i])," ",(evfvecs[i+1]+evfvecs[i])/2," ",1.0/(evfvecs[i+1]-evfvecs[i]))
       else
          println(ioa,(eigvs[1][i+2]+eigvs[1][i])/2," ",2.0/(eigvs[1][i+2]-eigvs[1][i])," ",(evfvecs[i+2]+evfvecs[i])/2," ",1.0/(evfvecs[i+2]-evfvecs[i])) 
       end    
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
   open("output/quasienergies_spacing.dat","w") do io
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


function buildingstate(lista,nn)
bc=FockBasis(nn)
ba=SpinBasis(1//2)
psif= 0.0*fockstate(bc, 0) ⊗ spindown(ba)
for i in 0:(nn-1)
    psidown = fockstate(bc, i) ⊗ spindown(ba)
    psiup   = fockstate(bc, i) ⊗ spinup(ba)
    psif = psif + lista[2*i+2]*psiup + lista[2*i+2-1]*psidown 
end
return psif
end

function buildingrho(rhoi,Nmax)
rho_ev = eigvals(rhoi)
rho_evecs = eigvecs(rhoi)
bc=FockBasis(Nmax)
ba=SpinBasis(1//2)
psi0= 0.0*fockstate(bc, 0) ⊗ spindown(ba)
rho = dm(psi0)  
for i in 1:length(rho_ev)
   listrho = [rho_evecs[j,i] for j in 1:length(rho_ev)]
   psiinst = buildingstate(listrho,Nmax)
   rhoinst = dm(psiinst)
   #println(rhoinst)
   rho = rho + rho_ev[i]*rhoinst
   #println(rho)
end
return rho
end


function negativitycorr(rho1,rho2,Nmax,L,r)
    xpint=0.05
    x = [-L:xpint:L;]
    rho1_qo = buildingrho(rho1,Nmax)
    rho1_pt = ptrace(rho1_qo,2)
    rho2_qo = buildingrho(rho2,Nmax)
    rho2_pt = ptrace(rho2_qo,2)
    sumwc   = 0.0
    sumwn1  = 0.0
    sumwn2  = 0.0
    sumwncor = 0.0
    sumwncorl = 0.0
    sumwncorl = 0.0
    sumIp  = 0.0
    sumIn  = 0.0
    sumwnd   = 0.0
    sumwndl  = 0.0
    sumwp1  = 0.0
    sumwp2  = 0.0
    sumwpcor = 0.0
    for ix in x
        for ip in x
           w1  = wigner(rho1_pt, ix, ip)
           w2 = wigner(rho2_pt, ix, ip)
           sumwc =  sumwc + xpint^2*real(w1)*real(w2)
           if real(w1)<0
               sumwn1 = sumwn1 + xpint^2*real(w1)^2
           else
               sumwp1 = sumwp1 + xpint^2*real(w1)^2
           end
           if real(w2)<0
               sumwn2 = sumwn2 + xpint^2*real(w2)^2
           else
               sumwp2 = sumwp2 + xpint^2*real(w2)^2
           end
           if real(w1) < -10^(-3) && real(w2) < -10^(-3)
               sumwncor = sumwncor + xpint^2*(real(w1)*real(w2))
               sumwnd   = sumwnd  + xpint^2*abs2(real(w1)-real(w2))
               if (ix^2/r+ip^2/r)<(1/r)
                   sumwncorl = sumwncorl + xpint^2*(real(w1)*real(w2))
                   sumwndl   = sumwndl   + xpint^2*abs2(real(w1)-real(w2))
               end    
           end
           if real(w1) > 0 && real(w2) > 0.0
               sumwpcor = sumwpcor + xpint^2*(real(w1)*real(w2))
           end
           if real(w1)*real(w2) < 0.0 
               sumIn = sumIn + xpint^2*real(w1)*real(w2)
           else
               sumIp = sumIp + xpint^2*real(w1)*real(w2)
           end    
        end
    end
    if  (sumwn1)>1e-8 && (sumwn2)>1e-8
        cn = sumwncor/(sumwn1*sumwn2)^(1/2)
    else
        cn = 10^(-10)
    end    
    cp = sumwpcor/(sumwp1*sumwp2)^(1/2)
    ov  = 2*pi*sumwc
    return [real(cn),real(cp),real(ov),sumwncor,sumwncorl,sumwnd,sumwndl,sumIp,sumIn]
end   

end
