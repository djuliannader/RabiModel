module DQPT
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
export amplitud


function initialstatequench(n,om,r,lambda,delta,eta,psi)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 evecs=eigvecs(HMatrix)
 gswf=[evecs[i,1] for i in 1:2*(n+1)]
 #println(gswf)
 return gswf
end



function amplitud(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi)
         #tint=2*pi/1.8
	 tint=0.1
	 tinti=0.01
	 nt=trunc(Int,tmax/tint)
	 nt2=40
	 t=0.0
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 eigvecs1=eigvecs(HMatrix)
	 eigvals1=eigvals(HMatrix)
	 #coef=[]
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 #for i in 1:length(psi0)
         #   suma=0.0
	 #   for j in 1:length(psi0)
	 #     suma = suma + conj(eigvecs1[j,i])*psi0[j]
	 #   end
         #   append!(coef,suma)
	 #end
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("Loschmidt_amplitud.dat","w") do io
 	 for i in 1:nt+1
 	     evol=exp(-im*HMatrix*t/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
	     #spp=0
	     #for k in 1:length(psi0)
	     #   spp = spp + exp(-im*eigvals1[k]*t)*abs2(coef[k]) 
             #end
	     #spp=abs2(spp)
 	     println(io,t," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
 	     t=t+tint
 	    end
 	 end
	 t=0.0
	 println("-------------   Go to file Loschmidt_amplitud.dat to see the the Loschmidt amplitud  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
	 open("Loschmidt_amplitud_ct.dat","w") do io
 	 for i in 1:nt+1
	   tim=-0.2
	   for j in 1:nt2+1
 	     evol=exp(-im*HMatrix*(t+tim*im)/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
	     #spp=0
	     #for k in 1:length(psi0)
	     #   spp = spp + exp(-im*eigvals1[k]*t)*abs2(coef[k]) 
             #end
	     #spp=abs2(spp)
 	     println(io,t," ",tim," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
	     tim=tim+tinti	     
 	    end
	    t=t+tint	     
	 end   
 	 end
	 println("-------------   Go to file Loschmidt_amplitud_ct.dat to see the the complex time Loschmidt amplitud  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
	 return "Done"
	 end

n=50
om=1.0
r=50.0
lambda0=0.0
delta0=0.0
eta0=0.0
psi0=0.0
lambda1=0.0
delta1=0.0
eta1=5.0
psi1=0.0


istate = initialstatequench(n,om,r,lambda0,delta0,eta0,psi0)
amplitudtest = amplitud(istate,10.0,1.0,n,om,r,lambda1,delta1,eta1,psi1)



end
