module DQPT
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import diagonalization
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


function survivalprobabilityt(rho0,tmax,n,om,r,lambda,delta,eta,psi)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tmax)
 tint=0.05
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 tinst=0.0
 nt = floor(Int, tmax/tint)
 surprob = []
 open("Loschmidt_amplitud_thermal.dat","w") do io
 for i in 1:(nt+1)
   rhot = sol(tinst)
   fidinst = (tr((rhot^(1/2)*rho0*rhot^(1/2))^(1/2)))
   println(io,tinst," ", round(real(fidinst),digits=16)," ",round(imag(fidinst),digits=16))
   tinst=tinst+tint
 end
 end
 return "done"
end

function survivalprobabilityt_ct(rho0,tmax,n,om,r,lambda,delta,eta,psi)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 nr=100
 ni=40
 inttre=tmax[1]/nr
 inttim=2*tmax[2]/ni
 open("Loschmidt_amplitud_thermal_ct.dat","w") do io
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
 return "done"
end



 
n=100
om=1.0
r=20.0
lambda0=0.0
delta0=0.0
eta0=2.0
psi0=0.0
eta1=1.125
lambda1=6.0
delta1=0.0
psi1=0.0

beta=10.0
tmax=10.0


rhoistate = initialthermalstate(n,om,r,lambda0,delta0,eta0,psi0,beta)

println("Size of the Fock basis: ",n)
println("Beta:                   ",beta)
sp = survivalprobabilityt(rhoistate,tmax,n,om,r,lambda1,delta1,eta1,psi1)
println(sp)
sp_c = survivalprobabilityt_ct(rhoistate,[tmax,0.5],n,om,r,lambda1,delta1,eta1,psi1)
println(sp_c)
println("Real-time survival amplitud data in file Loschmidt_amplitud_thermal.dat")

#ftest = (tr((rhoistate^(1/2)*rhoistate*rhoistate^(1/2))^(1/2)))^2
#println(ftest)

 
end