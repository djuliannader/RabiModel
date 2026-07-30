module DQPT_thermal
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
using QuantumOptics
include("diagonalization.jl")
include("wigner_eig.jl")
using .diagonalization
using .wigner_eig
export amplitud


function survivalprobabilityt(rho0,Jp,Jq,tlist,n,om,r,lambda,delta,eta,psi,L,tint,acc)
 tmax = tlist[length(tlist)]
 rho0qo = wigner_eig.buildingrho(rho0,n)
 rho0pt_b = ptrace(rho0qo,2)
 rho0pt_q = ptrace(rho0qo,1)
 bc = FockBasis(n)
 aop = destroy(bc)
 adop = dagger(aop)       
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tmax)
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=0.01,abstol =  acc)
 tinst=0.0
 nt = floor(Int, tmax/tint)
 surprob = []
 println("Density matrix obtained from t=0 to t=",tmax)  
 open("output/survival_probability_thermal.dat","w") do io
 open("output/negativities_quench_thermal.dat","w") do io2
 open("output/quantities_quench_thermal.dat","w") do io3        
 for tinst in tlist
   println(tinst)      
   rhot = sol(tinst)  
   neg = wigner_eig.wigner_negativities2(n,rhot,rho0,L,r)
   sp = tr(rhot*rho0)
   #mat = rho0*rhot
   #ls = [mat[i,i] for i in 1:(2*(n+1))]
   #lecho   = sum(ls)
   rhoqo = wigner_eig.buildingrho(rhot,n)
   rhopt_b = ptrace(rhoqo,2)
   rhopt_q = ptrace(rhoqo,1)
   sp_b  = tr(rho0pt_b*rhopt_b)
   sp_q  = tr(rho0pt_q*rhopt_q)
   fid   = (fidelity(rho0qo,rhoqo))^2
   fid_b = (fidelity(rho0pt_b,rhopt_b))^2  
   fid_q = (fidelity(rho0pt_q,rhopt_q))^2
   pur0 = tr(rho0qo^2)
   purt = tr(rhoqo^2)
   pur0_b = tr(rho0pt_b^2)
   purt_b = tr(rhopt_b^2)
   pur0_q = tr(rho0pt_q^2)
   purt_q = tr(rhopt_q^2)
   sf =  real(sp) + (1-pur0)^(1/2)*(1-purt)^(1/2)
   sf_b =  real(sp_b) + (1-pur0_b)^(1/2)*(1-purt_b)^(1/2)
   sf_q =  real(sp_q) + (1-pur0_q)^(1/2)*(1-purt_q)^(1/2)
   nav = real(tr(adop*aop*rhopt_b))
   n2av = real(tr((adop*aop)^2*rhopt_b))
   varn = real(n2av - nav)  
   println(io,tinst," ", real(sp)," ",real(sp_b)," ",real(sp_q)," ",real(fid)," ",real(fid_b)," ",real(fid_q)," ",real(sf)," ",real(sf_b)," ",real(sf_q))
   println(io2,tinst," ",(join(neg, " ")))
   #println(io3,tinst," ",nav," ",varn)    
   println(io3,tinst," ",nav/(varn)^(1/2)," ",varn/nav)  
 end
 end
 end
 end
 println("-------------   Go to file output/survival_probability_thermal.dat to see the Survival probability  ---------")
 println("             The file contains the survival probability from 0 to ",tmax," in steps of ",tint," time units ")
 println("-------------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativities_quench_thermal.dat to see the negativities     -----------")
 println("             The file contains Negativities from 0 to ",tmax," in steps of ",tint," time units           ")
 println("-------------------------------------------------------------------------------------------------------- ")
 return "done"
end

function wignerrhot(rho0,Jp,Jq,tshot,n,om,r,lambda,delta,eta,psi,L)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tshot)
 tint=0.05
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 rhot = sol(tshot)
 wig = wigner_eig.wigner_rhot(rhot,L,r,n)
 return "done"
end

function rhot(rho0,Jp,Jq,tshot,n,om,r,lambda,delta,eta,psi,L)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tshot)
 tint=0.05
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 rhot = sol(tshot)
 return rhot
end

function photonP(rho0,Jp,Jq,tshot,n,om,r,lambda,delta,eta,psi,L)
    rho = rhot(rho0,Jp,Jq,tshot,n,om,r,lambda,delta,eta,psi,L)
    rhoqo = wigner_eig.buildingrho(rho,n)
    rhopt_b = ptrace(rhoqo,2)
    bc = FockBasis(n)
    open("output/photondistribution.dat","w") do io
        for in in 0:20
            fockstaten = fockstate(bc,in)
            rhon = dm(fockstaten)
            rhonn = tr(rhopt_b*rhon)
            println(io,in," ",real(rhonn))
        end    
    end    
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
