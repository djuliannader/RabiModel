module dynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
export initialcoherent
export survivalp


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



function survivalp(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64)
         tint=0.05
	 nt=trunc(Int,tmax/tint)
	 t=0.0
	 # diagonal matrix elements
         vdiag=[(i)*om+r*(1/2)*(-1)^j for i in 0:Nmax for j in -1:0]
         #println(vdiag)
         HMatrix=Array(Diagonal(vdiag))
         # non diagonal matrix elements
         for i in 1:Nmax
            HMatrix[2i+1,2i]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
            HMatrix[2i,2i+1]=lambda*((delta+1)/2)*(i)^(1/2)   # Jaynes-Cummings
            HMatrix[2+2i,2i-1]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
            HMatrix[2i-1,2+2i]=lambda*((1-delta)/2)*(i)^(1/2)   # Anti Jaynes-Cummings
         end
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0)) 
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("survivalprobability.dat","w") do io
 	 for i in 1:nt+1
 	     evol=exp(-im*HMatrix*t/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     sp=abs2(sp[1])
 	     println(io,t," ",round(sp,digits=16))
 	     t=t+tint
 	    end
 	 end
	 println("-------------   Go to file survivalprobability.dat to see the survival probability  ----------------")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
	 return "Done"
	 end

function survivalpt(psi0::Vector{Complex{Float64}},fq::Matrix{Complex{Float64}},tmax::Float64,om)
 pi=acos(-1)
 T=2*pi/om
 nint=trunc(Int64,tmax/T)
 psi0a=Array{Complex{Float64}}(undef,1,length(psi0)) 
 for k in 1:length(psi0)
     psi0a[1,k]=conj(psi0[k])
 end
 psi0t=psi0
 open("survivalprobability_f.dat","w") do io
 for i in 0:nint
   sp=psi0a*psi0t
   sp=abs2(sp[1])
   println(io,T*i," ",round(sp,digits=16))
   psi0t=fq*psi0t
 end
 end
 println("-------------   Go to file survivalprobability_f.dat to see the survival probability  --------------")
 println("             The file contains SP from 0 to ",tmax," in steps of ",T ," time units               ")
 println("--------------------------------------------------------------------------------------------------- ")
 return "done"
end


#test=initialcoherent(1.0,0.0,1.0,1.5,1.0,50)
#mensaje=survivalp(test,50.0,1.0,50,1.0,1.0,0.05,0.0)
#println("Prueba: ",mensaje)


end