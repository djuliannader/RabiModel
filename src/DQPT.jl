module DQPT
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
export amplitud
export initialstatequench
export overlapdqpt
export Nzeros
export PositionsZeros
export overlapdqptlist
import troterization


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
   open("overlap_dqpt.dat","w") do io
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



function amplitud(psi0::Vector{Complex{Float64}},tmax::Float64,hbar::Float64,Nmax::Int64,om::Float64,r::Float64,lambda::Float64,delta::Float64,eta,psi)
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
	 open("Loschmidt_amplitud.dat","w") do io
 	 for i in 1:nt+1
 	     evol=exp(-im*HMatrix*t/hbar)
	     psi0t=evol*psi0
 	     sp=psi0a*psi0t
 	     spf=sp[1]
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
	 println("-------------   Go to file Loschmidt_amplitud_ct.dat to see the the complex time Loschmidt amplitud  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
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
	 println("-------------   Go to file Loschmidt_amplitud_fld.dat to see the the Loschmidt amplitud  ----------------")
	 println("----------- Dynamics governed by the Floquet operator for the time dependent Hamiltonian       -----")
	 println("             The file contains SP for ",ntmax," period steps of ",tint," time units               ")
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
	 println("-------------   Go to file Loschmidt_amplitud_fld_ct.dat to see the the complex time Loschmidt amplitud  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",ntmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
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
	 println("-------------   Go to file Jz.dat to see the the time evolution of the Jz operator  ----------------")
	 println("----------- Dynamics governed by the time evolution operator for the time independent Hamiltonian       -----")
	 println("             The file contains SP from 0 to ",tmax," in steps of ",tint," time units               ")
	 println("--------------------------------------------------------------------------------------------------- ")
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


#n=110
#om=1.0
#r=20.0
#lambda0=0.0
#delta0=0.0
#eta0=3/2
#psi0=0.0
#eta1=3/2
#lambda1=1/5
#delta1=0.0
#psi1=0.0
#nn=500
#nsubint=1000
#nsubint2=400
#nsubint3=30

#name = "position_zeros.dat"

#alpha=1.0
#ph=0.0

#tmax=10.0
#tcirc=[0.0-0.5*im,0.0+0.5*im,10.0+0.5*im,10.0-0.5*im]
#tcircr=[0.0-0.0*im,0.0+0.5*im,10.0+0.5*im,10.0-0.0*im]
#tcircl=[0.0-0.5*im,0.0+0.0*im,10.0+0.0*im,10.0-0.5*im]

#istate = initialstatequench(n,om,r,lambda0,delta0,eta0,psi0)

#phi0 = alpha^(1/2)*istate[1] + (1-alpha)^(1/2)*exp(im*ph)*istate[2]

#hamf = diagonalization.hamiltonian(n,om,r,lambda1,delta1,eta1,psi1)

#amplitudint = amplitud(phi0,tmax,1.0,n,om,r,lambda1,delta1,eta1,psi1)






#computingjz = Jz(phi0,tmax,1.0,n,om,r,lambda1,delta1,eta1,psi1)
#ovl = overlapdqpt(phi0,hamf)
#println("Overlap Done")
#rr = Nzeros(phi0,tcirc,1.0,n,om,r,lambda1,delta1,eta1,psi1,nsubint)
#rrr = Nzeros(phi0,tcircr,1.0,n,om,r,lambda1,delta1,eta1,psi1,nsubint)
#rrl = Nzeros(phi0,tcircl,1.0,n,om,r,lambda1,delta1,eta1,psi1,nsubint)
#println("Number of zeros within the full circuit: ",rr)
#println("Number of zeros on the right of the real axis: ",rrr)
#println("Number of zeros on the left  of the real axis: ",rrl)
#println("---------------------------------------------------------")

#println("-Calculating the position of the zeros in the complex plane-")
#pos = PositionsZeros(phi0,tcirc,1.0,n,om,r,lambda1,delta1,eta1,psi1,nsubint2,nsubint3,name)

#println("Number of zeros found in the contour rectangle: ",pos)
#println("Their positions in the complex plane can be found in file position_zeros.dat")

end
