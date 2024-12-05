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
 println("here")
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



 
n=50
om=1.0
r=10.0
lambda0=0.5
delta0=0.0
eta0=2.0*r^(1/2)
psi0=0.0
eta1=(2/3)*r^(1/2)
lambda1=0.5
delta1=0.0
psi1=0.0

beta=2.0
tmax=10.0


rhoistate = initialthermalstate(n,om,r,lambda0,delta0,eta0,psi0,beta)

sp = survivalprobabilityt(rhoistate,tmax,n,om,r,lambda1,delta1,eta1,psi1)
 
#ftest = (tr((rhoistate^(1/2)*rhoistate*rhoistate^(1/2))^(1/2)))^2
#println(ftest)

 
end