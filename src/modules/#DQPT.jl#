module DQPT
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("diagonalization.jl")
include("wigner_eig.jl")
include("Fisher.jl")
include("troterization.jl")
using .diagonalization
using .Fisher
using .wigner_eig
using .troterization
export amplitud
export initialstatequench
export overlapdqpt
export Nzeros
export PositionsZeros
export overlapdqptlist



function initialstatequench(n,om,r,lambda,delta,eta,psi)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 evecs=eigvecs(HMatrix)
 gswf=[evecs[i,1] for i in 1:2*(n+1)]
 eswf=[evecs[i,2] for i in 1:2*(n+1)]
 #println(gswf)
 return [gswf,eswf]
end


function overlapdqpt(psi0,Hamf,nmax)
   eners =eigvals(Hamf) 
   states=eigvecs(Hamf)
   open("output/overlap_dqpt.dat","w") do io
   for i in 1:length(eners)
     istate=[states[j,i] for j in 1:2*(nmax+1)]
     psi0t = transpose(conj(psi0))
     overlap2 = abs2(psi0t*istate)
     println(io,eners[i]," ",overlap2)
   end
   end
end

function overlapdqptlist(psi0,Hamf,nmax)
   listovlp=[]
   states=eigvecs(Hamf)
   for i in 1:length(psi0)
     istate=[states[j,i] for j in 1:2*(nmax+1)]
     psi0t = transpose(conj(psi0))
     overlap2 = abs2(psi0t*istate)
     append!(listovlp,overlap2)
   end
   return listovlp
end



function amplitud(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi,L)
	 tint=0.05
	 tinti=0.01
	 nt=trunc(Int,tmax/tint)
	 nt2=100
	 t=0.0
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 eigvecs1=eigvecs(HMatrix)
	 eigvals1=eigvals(HMatrix)
	 #println("ev1:",eigvals1[1])
	 # adjunt initial state
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("output/Loschmidt_amplitud.dat","w") do io
	 open("output/negativities_quench.dat","w") do io2
	 open("output/qfi_quench.dat","w") do io3
	 open("output/zeros_cutwigner.dat","w") do io4
 	 for i in 1:nt+1
 	     evol=exp(-im*HMatrix*t/hbar)
	     psi0t=evol*psi0
	     neg = wigner_eig.wigner_negativities(Nmax,psi0t,L)
	     rho =  psi0t*transpose(conj(psi0t))
	     zeros =  wigner_eig.wigner_rhot_ZerosCut(rho,Nmax,r)
	     rhoqo = wigner_eig.buildingrho(rho,Nmax)
	     rhopt = ptrace(rhoqo,2)
	     qfi1 = Fisher.fishern2(rhopt,Nmax)
	     qfi2 = Fisher.fisherdisplacementx(rhopt,Nmax)
             qfi3 = Fisher.fisherdisplacementp(rhopt,Nmax)
 	     sp=psi0a*psi0t
 	     spf=sp[1]
 	     println(io,t," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
	     println(io2,t," ",round(real(neg),digits=8))
	     println(io3,t," ",round(real(qfi1),digits=8)," ",round(real(qfi2),digits=8)," ",round(real(qfi3),digits=8))
	     println(io4,t," ",round(real(zeros),digits=8))
 	     t=t+tint
 	    end
	 end
	 end
	 end
	 end
	 t=0.0
	 println("------------------------------------------------------------------------------------------------------------ ")
	 println("-------------   Go to file output/Loschmidt_amplitud.dat to see the Loschmidt amplitud as function of time---")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains the Survival amplitude from 0 to ",tmax," in steps of ",tint," time units     ")
	 println("------------------------------------------------------------------------------------------------------------ ")
	 println("-------------   Go to file output/negativities_quench.dat to see the negativity volume as a function of time-")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains the negativity volume from t=0 to ",tmax," in steps of ",tint," time units    ")
	 println("------------------------------------------------------------------------------------------------------------ ")
	 println("-------------   Go to file output/qfi_quench.dat to see the QFI as a function of time   ---------------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains the negativity volume from t=0 to ",tmax," in steps of ",tint," time units    ")
	 println("------------------------------------------------------------------------------------------------------------ ")
	 println("-------------   Go to file output/zeros_cutwigner.dat to the the number of sign changes (Ns)       ----------")
	 println("-------------   of the Wigner function cuts along the position and momentum axis as a function of time ------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains Ns from 0 to ",tmax," in steps of ",tint," time units                         ")
	 println("------------------------------------------------------------------------------------------------------------ ")
 	 open("output/Loschmidt_amplitud_ct.dat","w") do io
 	 for i in 1:nt+1
	   tim=-0.5
	   for j in 1:nt2+1
 	     evol=exp(-im*HMatrix*(t+tim*im)/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
 	     println(io,t," ",tim," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
	     tim=tim+tinti	     
 	    end
	    t=t+tint	     
	 end   
 	 end
	 println("-------------   Go to file output/Loschmidt_amplitud_ct.dat to see the  complex time Loschmidt amplitud  -----")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("------------------------------------------------------------------------------------------------------------ ")
	 return "Done"
	 end


function amplitud2(psi0::Vector{Complex{Float64}},ntmax::Integer,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,gamma,omega,eta,psi,nn)
         #tint=(2*pi/gamma)*(1/100)
	 tint=(2*pi/gamma)*(1)	 
	 # building the evolution operator
         floquet = troterization.troter2(Nmax,nn,r,om,gamma,omega,eta,psi)
	 floquet_ct = troterization.troter2_ct(Nmax,nn,r,om,gamma,omega,eta,psi)
	 #floquetm= floquet^(1/100)
	 #floquetm2= floquet_ct^(1/100)
	 # conjugate initial state
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("Loschmidt_amplitud_fld.dat","w") do io
	 println(io,0.0," ",1.0," ", 0.0)
 	 for i in 1:ntmax
	     psi0t=floquet*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
	     t=tint*i
 	     println(io,t," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
	     psi0=psi0t
 	 end
	 end
	 println("-------------   Go to file output/Loschmidt_amplitud_fld.dat to see the the Loschmidt amplitud  ----")
	 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
	 println("             The file contains SP for ",ntmax," period steps of ",tint," time units                 ")
	 println("--------------------------------------------------------------------------------------------------- ")
	 open("Loschmidt_amplitud_fld_ct.dat","w") do io
         for i in 0:ntmax 
           for j in -5:5
	     floquet1p = floquet^(i)
	     floquet2p = floquet_ct^(j)
	     psi0t=floquet2p*floquet1p*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
	     t=tint*i
	     tint*j
 	     println(io,tint*i," ",tint*j," ",round(real(spf),digits=16)," ", round(imag(spf),digits=16))
 	   end
	 end
	 end
	 println("-------------   Go to file output/Loschmidt_amplitud_fld_ct.dat to see the the complex time Loschmidt amplitud  ----")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       ------------")
	 println("             The file contains SP from 0 to ",ntmax," in steps of ",tint," time units                               ")
	 println("------------------------------------------------------------------------------------------------------------------- ")
	 return "Done"
	 end


     function Jz(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi)
	 tint=0.02
	 nt=trunc(Int,tmax/tint)
	 t=0.0
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 diagJz=[r*(1/2)*(-1)^j+0*im for i in 0:n for j in -1:0]
	 Jz=Array(Diagonal(diagJz))
	 eigvecs1=eigvecs(HMatrix)
	 eigvals1=eigvals(HMatrix)
	 #println("ev1:",eigvals1[1])
	 # adjunt initial state
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 open("Jz.dat","w") do io
 	 for i in 1:nt+1
 	     jzevol=exp(im*HMatrix*t/hbar)*Jz*exp(-im*HMatrix*t/hbar)
	     psi0jt=jzevol*psi0
 	     jzt=psi0a*psi0jt
 	     jzev=jzt[1]+r/2
 	     println(io,t," ",round(real(jzev),digits=16))
 	     t=t+tint
 	    end
 	 end
	 t=0.0
	 println("-------------   Go to file output/Jz.dat to see the the time evolution of the Jz operator  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       ---")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("---------------------------------------------------------------------------------------------------------- ")
	 return "Done"
	 end



function Nzeros(psi0::Vector{Complex{Float64}},tcirc::Vector{Complex{Float64}},hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi,nsubint)
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 evf=eigvecs(HMatrix)
	 evals = eigvals(HMatrix)
         # adjunt initial state
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 coef2=[]
	 for k in 1:length(psi0)
	   listevf= [evf[i,k] for i in 1:length(psi0)]
	   ck= psi0a*listevf
	   append!(coef2,abs2(ck[1]))
	 end
	 # extracting intervals
	 nint = nsubint
         tint1 = abs((imag(tcirc[2])-imag(tcirc[1]))/nint)
	 tint2 = abs((real(tcirc[3])-real(tcirc[2]))/nint)
	 tint3 = abs((imag(tcirc[4])-imag(tcirc[3]))/nint)
	 tint4 = abs((real(tcirc[1])-real(tcirc[4]))/nint)
	 tinst = tcirc[1]
	 println("#--- corners of the rectangle for the contour integral ----#")
	 println("(0) t: ",tinst)
#	 segment 1 (a lo largo del eje imaginario)
         sum1 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum1 = sum1 + im*fz*tint1 
            tinst = tinst + tint1*im
	 end
	 println("(1) t:", tinst)
#	 segment 2 (a lo largo del eje real)
         sum2 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum2 = sum2 + fz*tint2 
            tinst = tinst + tint2
	 end
	 println("(2) t: ",tinst)
#	 segment 3 (a lo largo del eje imaginario)
         sum3 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum3 = sum3 - fz*tint3*im 
            tinst = tinst - tint3*im
	 end
	 println("(3) t: ",tinst)
#	 segment 4 (a lo largo del eje real)
         sum4 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum4 =  sum4 - fz*tint4 
            tinst = tinst - tint4
	 end
	 println("(4) t: ",tinst)
         println("#-----------------------------------------------#")
	 sumt=imag(sum1+sum2+sum3+sum4)
         return -1.0*sumt/(2*pi) 
	 end

function Nzeros2(psi0::Vector{Complex{Float64}},tcirc::Vector{Complex{Float64}},hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi,nsubint)
	 # building Hamiltonian
         HMatrix= diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
	 evf=eigvecs(HMatrix)
	 evals = eigvals(HMatrix)
         # adjunt initial state
	 psi0a=Array{Complex{Float64}}(undef,1,length(psi0))
	 for k in 1:length(psi0)
	     psi0a[1,k]=conj(psi0[k])
         end
	 coef2=[]
	 for k in 1:length(psi0)
	   listevf= [evf[i,k] for i in 1:length(psi0)]
	   ck= psi0a*listevf
	   append!(coef2,abs2(ck[1]))
	 end
	 # extracting intervals
	 nint = nsubint
         tint1 = abs((imag(tcirc[2])-imag(tcirc[1]))/nint)
	 tint2 = abs((real(tcirc[3])-real(tcirc[2]))/nint)
	 tint3 = abs((imag(tcirc[4])-imag(tcirc[3]))/nint)
	 tint4 = abs((real(tcirc[1])-real(tcirc[4]))/nint)
	 tinst = tcirc[1]
#	 println("#--- corners of the rectangle for the contour integral ----#")
#	 println("(0) t: ",tinst)
#	 segment 1 (a lo largo del eje imaginario)
         sum1 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum1 = sum1 + im*fz*tint1 
            tinst = tinst + tint1*im
	 end
#	 println("(1) t:", tinst)
#	 segment 2 (a lo largo del eje real)
         sum2 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum2 = sum2 + fz*tint2 
            tinst = tinst + tint2
	 end
#	 println("(2) t: ",tinst)
#	 segment 3 (a lo largo del eje imaginario)
         sum3 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum3 = sum3 - fz*tint3*im 
            tinst = tinst - tint3*im
	 end
#	 println("(3) t: ",tinst)
#	 segment 4 (a lo largo del eje real)
         sum4 = 0.0 + 0.0*im
         for it in 1:nint
            list1 = [coef2[k]*evals[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
	    list2 = [coef2[k]*exp(-im*tinst*evals[k]) for k in 1:length(psi0)]
            fz = (1/sum(list2))*(-im)*sum(list1)
	    sum4 =  sum4 - fz*tint4 
            tinst = tinst - tint4
	 end
#	 println("(4) t: ",tinst)
#         println("#-----------------------------------------------#")
	 sumt=imag(sum1+sum2+sum3+sum4)
         return -1.0*sumt/(2*pi) 
	 end


function PositionsZeros(psi0::Vector{Complex{Float64}},tcirc::Vector{Complex{Float64}},hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi,nsubint,nsubint2,name)
         deltat1 = abs((imag(tcirc[2])-imag(tcirc[1]))/nsubint)
	 deltat2 = abs((real(tcirc[3])-real(tcirc[2]))/nsubint)
	 tcircint = [tcirc[1],tcirc[1]+deltat1*im,tcirc[1]+deltat1*im +deltat2, tcirc[1]+deltat2]
	 counter = 0
	 open(name,"w") do io
	 for i in 1:(nsubint)
	   for j in 1:(nsubint)
	      zero = Nzeros2(psi0,tcircint,1.0,Nmax,om,r,lambda,delta,eta,psi,nsubint2)
	      if round(Int64, zero)>0
	       counter = counter + 1
	       println(counter," ",(imag(tcircint[1])+imag(tcircint[2]))/2.0," ",(real(tcircint[1])+real(tcircint[3]))/2.0)
	       println(io,(imag(tcircint[1])+imag(tcircint[2]))/2.0," ",(real(tcircint[1])+real(tcircint[3]))/2.0)
	      end
	      tcircint[1]=real(tcircint[1])+im*(imag(tcircint[1])+deltat1)
	      tcircint[2]=real(tcircint[2])+im*(imag(tcircint[2])+deltat1)
	      tcircint[3]=real(tcircint[3])+im*(imag(tcircint[3])+deltat1)
	      tcircint[4]=real(tcircint[4])+im*(imag(tcircint[4])+deltat1)
	   end
	   tcircint[1]=real(tcirc[1])+i*deltat2+im*(imag(tcirc[1]))
	   tcircint[2]=real(tcirc[1])+i*deltat2+im*(imag(tcirc[1])+deltat1)
	   tcircint[3]=real(tcirc[1])+(i+1)*deltat2+im*(imag(tcirc[1])+deltat1)
	   tcircint[4]=real(tcirc[1])+(i+1)*deltat2+im*(imag(tcirc[1]))
	 end
	 end
         return counter
end


function wigner_rhot(psi0,ham,L,r,Nmax,t)
 psit = exp(-im*t*ham)*psi0
 rhot = psit*transpose(conj(psit))
 wig = wigner_eig.wigner_rhot(rhot,L,r,Nmax)
 return "done"
end



end
