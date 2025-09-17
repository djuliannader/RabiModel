module DQPT_thermal
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("diagonalization.jl")
include("wigner_eig.jl")
using .diagonalization
using .wigner_eig
export amplitud


function initialthermalstate(n,om,r,lambda,delta,eta,psi,beta)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 evecs=eigvecs(HMatrix)
 evals=eigvals(HMatrix)
 rho0= zeros(Complex, length(evals), length(evals))
 facs=[]
 for nn in 1:length(evals)
   fac=exp(-(evals[nn]-1.0*evals[1])*beta)
   append!(facs,fac)
 end
 facsn=(1/sum(facs))*facs
 for nn in 1:length(evals)
   vecn   = [evecs[i,nn] for i in 1:length(evals)]
   vecntp = transpose(conj(vecn))
   rho0 = rho0 + facsn[nn]*vecn*vecntp
 end
 diagnumop=[(i)+0.0*im for i in 0:n for j in -1:0]
 numop=Diagonal(diagnumop)
 expnum = tr(numop*rho0)
 println("Initial expectation value of the number operator: ",expnum)
 return rho0
end


function survivalprobabilityt(rho0,tmax,n,om,r,lambda,delta,eta,psi,L)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tmax)
 tint=0.01
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 tinst=0.0
 nt = floor(Int, tmax/tint)
 surprob = []
 open("output/Loschmidt_amplitud_thermal.dat","w") do io
 open("output/negativities_quench_thermal.dat","w") do io2
 for i in 1:(nt+1)
   rhot = sol(tinst)
   neg = wigner_eig.wigner_rhot_neg(rhot,n,L)
   lecho = (tr((rhot^(1/2)*rho0*rhot^(1/2))^(1/2)))
   #mat = rho0*rhot
   #ls = [mat[i,i] for i in 1:(2*(n+1))]
   #lecho   = sum(ls)
   println(io,tinst," ", round(real(lecho),digits=8)," ",round(imag(lecho),digits=8))
   println(io2,tinst," ", round(real(neg),digits=8))
   tinst=tinst+tint
 end
 end
 end
 println("-------------   Go to file output/Loschmidt_amplitud_thermal.dat to see the Loschmidt amplitud  ---------")
 println("             The file contains the survival amplitude from 0 to ",tmax," in steps of ",tint," time units ")
 println("-------------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativities_quench_thermal.dat to see the negativities     -----------")
 println("             The file contains Negativities from 0 to ",tmax," in steps of ",tint," time units           ")
 println("-------------------------------------------------------------------------------------------------------- ")
 return "done"
end

function wignerrhot(rho0,tshot,n,om,r,lambda,delta,eta,psi,L)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tshot)
 tint=0.05
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 rhot = sol(tshot)
 wig = wigner_eig.wigner_rhot(rhot,L,r,n)
 return "done"
end

function survivalprobabilityt_ct(rho0,tmax,n,om,r,lambda,delta,eta,psi)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 nr=100
 ni=40
 inttre=tmax[1]/nr
 inttim=2*tmax[2]/ni
 open("output/Loschmidt_amplitud_thermal_ct.dat","w") do io
 tre=0
 for i in 1:(nr+1)
   tim=-tmax[2]
   for j in 1:(ni+1)
     rhot=exp(-im*HMatrix*(tre+im*tim))*rho0*exp(im*HMatrix*(tre-im*tim))
     #rhot=exp(im*HMatrix*(tre))*rho0*exp(-im*HMatrix*(tre))
     fidinst = (tr((rho0^(1/2)*rhot*rho0^(1/2))^(1/2)))
     println(io,tre," ",tim," ", round(real(fidinst),digits=16)," ",round(imag(fidinst),digits=16))
     #println(io,tre," ", round(real(fidinst),digits=16)," ",round(imag(fidinst),digits=16))
     tim=tim+inttim
   end
   tre=tre+inttre
 end
 end
 println("-------------   Go to file output/Loschmidt_amplitud_thermal.dat to see the complex-time Loschmidt amplitud  ---------")
 println("             The file contains the survival amplitude from 0 to ",tmax," in steps of ",tint," time units              ")
 println("--------------------------------------------------------------------------------------------------------------------- ")
 return "done"
end
 


 
end
