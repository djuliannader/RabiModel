module dynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
import diagonalization
import wigner_eig
import Fisher
export initialcoherent
export survivalp


function initialcat(xi::Float64,pii::Float64,ti::Float64,fi::Float64,hbar::Float64,Nmax::Int64)
         cat1=initialcoherent(xi,pii,ti,fi,hbar,Nmax)
         cat2=initialcoherent(-xi,pii,ti,fi,hbar,Nmax)
         catnm=cat1+cat2
	 catnmt=transpose(catnm)
	 nfac=catnmt*catnm
	 cat=(1/nfac^(1/2))*catnm
	 return cat
	 end

function initialsqueezed(xi::Float64,pii::Float64,ti::Float64,fi::Float64,hbar::Float64,Nmax::Int64,r::Float64,theta::Float64)
         cs = initialcoherent(xi,pii,ti,fi,hbar,Nmax)
	 aop=zeros(ComplexF64,2*(Nmax+1), 2*(Nmax+1))
	 for i in 1:Nmax
	    aop[(2*i-1),(2i+1)]=sqrt(i)
	    aop[2*i,2i+2]=sqrt(i)
	 end
	 adop=transpose(aop)
	 expsq = (r/2)*((cos(2*theta)-im*sin(2*theta))*aop^2-(cos(2*theta)+im*sin(2*theta))*adop^2) 
	 sqzop = exp(expsq)
	 sqzstate = sqzop*cs
	 #println(sqzstate)
	 sqzstatet=transpose(sqzstate)
	 nfac= sqrt(sqzstatet*sqzstate)
	 nsqzstate = (1/nfac)*sqzstate
	 #fock1=adop*cs
	 return nsqzstate
	 end


function initialcoherent(xi::Float64,pii::Float64,ti::Float64,fi::Float64,hbar::Float64,Nmax::Int64)
	 alabs=((xi^2+pii^2)*(1/(2*hbar)))^(1/2)
	 al=(1/(2*hbar)^(1/2))*(xi+im*pii)
	 z=tan(ti/2.0)*exp(-im*fi)
	 zabs2=abs2(z)
	 gcs=[(al^n)/((factorial(big(n)))^(1/2)) for n in 0:Nmax]
	 bcs=[((binomial(1,n))^(1/2))*z^(n) for n in 0:1]
	 gcs = exp(-alabs^2/2)*gcs
	 bcs = (1/(1+zabs2)^(1/2))*bcs
	 gcs=[round(gcs[i],digits=20) for i in 1:length(gcs)]
	 gcsf=[convert(Complex{Float64},gcs[i]) for i in 1:length(gcs)]
	 csf=[gcsf[i+1]*bcs[j] for i in 0:Nmax for j in 1:2]
	 return csf
	 end



function survivalp(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi)
         #tint=2*pi/1.8
	 tint=0.05
	 nt=trunc(Int,tmax/tint)
	 t=0.0
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0)) 
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("output/survivalprobability.dat","w") do io
 	 for i in 1:nt+1
 	     evol=exp(-im*HMatrix*t/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     sp=abs(sp[1])
 	     println(io,t," ",round(sp,digits=16))
 	     t=t+tint
 	    end
 	 end
	 println("-------------   Go to file output/survivalprobability.dat to see the survival probability  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
	 return "Done"
	 end

function survivalpt(psi0::Vector{Complex{Float64}},fq::Matrix{Complex{Float64}},tmax::Float64,om)
 pi=acos(-1)
 T=(2*pi/om)/1
 nint=trunc(Int64,tmax/T)
 psi0a=Array{Complex{Float64}}(undef,1,length(psi0)) 
 for k in 1:length(psi0)
     psi0a[1,k]=conj(psi0[k])
 end
 psi0t=psi0
 open("output/survivalprobability_f.dat","w") do io
 for i in 0:nint
   sp=psi0a*psi0t
   sp=abs(sp[1])
   println(io,T*i," ",round(sp,digits=16))
   psi0t=(fq^(1))*psi0t
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/survivalprobability_f.dat to see the survival probability  --------------")
 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
 println("             The file contains SP from 0 to ",tmax," in steps of ",T ," time units               ")
 println("--------------------------------------------------------------------------------------------------- ")
 return "done"
end

function fotoc(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi,L)
 tint=0.05
 t=0
 qfilist=[]
 neglist=[]
 nt=trunc(Int,tmax/tint)
 HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
 bc=FockBasis(Nmax)
 adop=create(bc)
 aop = destroy(bc)
 xop=(1/(2)^(1/2))*(aop+adop)
 pop=(im/(2)^(1/2))*(adop-aop)
 pi=acos(-1)
 T=(2*pi/om)/1
 nint=trunc(Int64,tmax/T)
 psit=psi0
 open("output/fotoc.dat","w") do io
 open("output/qfi.dat","w") do io2
 open("output/negativities.dat","w") do io3
 for i in 1:nt+1
   psit=exp(-im*HMatrix*t)*psi0
   neg = wigner_eig.wigner_negativities(Nmax,psit,L)
   phi = wigner_eig.buildingstate(psit,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   qfi=Fisher.fishern(rhopt,Nmax)
   x1m=expect(xop,rhopt)
   x2m=expect(xop^2,rhopt)
   p1m=expect(pop,rhopt)
   p2m=expect(pop^2,rhopt)
   fotoc = (x2m-x1m^2)^(1/2) + (p2m-p1m^2)^(1/2)
   println(io,t," ",round(real(fotoc),digits=16))
   println(io2,t," ",round(real(qfi),digits=8))
   println(io3,t," ",round(neg,digits=8))
   append!(qfilist,round(real(qfi),digits=8))
   append!(neglist,round(real(neg),digits=8))
   t=t+tint
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fotoc.dat to see the fotoc                           -------------")
 println("----------- Dynamics governed by the time independent Hamiltonian                              -----")
 println("             The file contains Fotoc(t) from t=0 to ",tmax," in steps of ",tint," time units        ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi.dat to see the Quantum Fisher Information       --------------")
 println("----------- Dynamics governed by the time independent Hamiltonian                              -----")
 println("             The file contains QFI(t) from t=0 to ",tmax," in steps of ",tint," time units          ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativities.dat to see the Wigner negativities       ------------")
 println("----------- Dynamics governed by the time independent Hamiltonian                              -----")
 println("             The file contains Neg(t) from t=0 to ",tmax," in steps of ",tint," time units          ")
 println("--------------------------------------------------------------------------------------------------- ")
 qfiavlist=[tint*qfilist[i] for i in 3:length(qfilist)]
 negavlist=[tint*neglist[i] for i in 1:length(neglist)]
 avqfi=(1/(tmax-2*tint)*sum(qfiavlist))
 avneg=(1/(tmax)*sum(negavlist))
 return [avqfi,avneg]
end



function fotoct(psi0::Vector{Complex{Float64}},fq::Matrix{Complex{Float64}},tmax::Float64,om::Float64,Nmax,L)
 bc=FockBasis(Nmax)
 adop=create(bc)
 aop = destroy(bc)
 xop=(1/(2)^(1/2))*(aop+adop)
 pop=(im/(2)^(1/2))*(adop-aop)
 pi=acos(-1)
 T=(2*pi/om)/1
 nint=trunc(Int64,tmax/T)
 psi0t=psi0
 qfilist=[]
 neglist=[]
 open("output/fotoc_f.dat","w") do io
 open("output/qfi_f.dat","w") do io2
 open("output/negativities_f.dat","w") do io3
 for i in 0:nint
   phi = wigner_eig.buildingstate(psi0t,Nmax)
   neg = wigner_eig.wigner_negativities(Nmax,psi0t,L)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   qfi=Fisher.fishern(rhopt,Nmax)
   x1m=expect(xop,rhopt)
   x2m=expect(xop^2,rhopt)
   p1m=expect(pop,rhopt)
   p2m=expect(pop^2,rhopt)
   fotoc = (x2m-x1m^2)^(1/2) + (p2m-p1m^2)^(1/2)
   println(io,T*i," ",round(real(fotoc),digits=8))
   println(io2,T*i," ",round(real(qfi),digits=8))
   println(io3,T*i," ",round(neg,digits=8))
   append!(qfilist,round(real(qfi),digits=8))
   append!(neglist,round(real(neg),digits=8))
   psi0t=(fq^(1))*psi0t
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fotoc_f.dat to see the fotoc                        --------------")
 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
 println("             The file contains SP from 0 to ",tmax," in steps of ",T ," time units               ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_f.dat to see the Quantum Fisher Information            -------")
 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
 println("             The file contains QFI(t) from t=0 to ",tmax," in steps of ",T ," time units          --")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativities_f.dat to see the Wigner negativities        ---------")
 println("----------- Dynamics governed by the the Floquet operator for the time dependent Hamiltonian   -----")
 println("             The file contains Neg(t) from t=0 to ",tmax," in steps of ",T ," time units          --")
 println("--------------------------------------------------------------------------------------------------- ")
 qfiavlist=[qfilist[i] for i in 3:length(qfilist)]
 negavlist=[neglist[i] for i in 1:length(neglist)]
 avqfi=(1/(length(qfiavlist))*sum(qfiavlist))
 avneg=(1/(length(negavlist))*sum(negavlist))
 return [avqfi,avneg]
end

function survivalpt2(psi0::Vector{Complex{Float64}},fq::Matrix{Complex{Float64}},tmax::Float64,om)
 pi=acos(-1)
 T=2*pi/om
 nint=trunc(Int64,tmax/T)
 psi0a=Array{Complex{Float64}}(undef,1,length(psi0)) 
 for k in 1:length(psi0)
     psi0a[1,k]=conj(psi0[k])
 end
 psi0t=psi0
 open("output/survivalprobability_f2.dat","w") do io
 for i in 0:nint
   sp=psi0a*psi0t
   sp=abs(sp[1])
   println(io,T*i," ",round(sp,digits=16))
   psi0t=fq*psi0t
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/survivalprobability_f2.dat to see the survival probability  ------")
 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
 println("             The file contains SP from 0 to ",tmax," in steps of ",T ," time units             -----")
 println("--------------------------------------------------------------------------------------------------- ")
 return "done"
end




end