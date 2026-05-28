module dqpts_thermal_tprime
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
using DifferentialEquations
include("diagonalization.jl")
include("DQPT_thermal.jl")
include("DQPT.jl")
include("wigner_eig.jl")
using .diagonalization
using .DQPT_thermal
using .DQPT
using .wigner_eig



n=70              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=10.0            # Qubit frequency
lambda0=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g0=0              # Initial coxupling
phs=0.0           # Phase of the Hamiltonian
g1=5/4            # Final coupling
lambda1=0.0       # Final Carrier parameter
tmax=10.0         # maximal time for the survival probability
tprime=0.0        # time for Wigner function
flag1=0           # (1) for the complex time survival probability (0) for skip 
beta=10.0         # temperature
L=7.5            # Size of the phase space
tint = 0.05       # time intervals
acc = 10^(-10)    # accuracy for the diferential equation
kp = 0.0          # photon loss
kq = 0.0          # qubit relaxation    



#------------- Perform calculations---------------------------#

println("Size of the Fock basis: ",n)
println("Beta:                   ",beta)


#----- Building initial state
# initial thermal state
if beta<100
    rhoistate = DQPT_thermal.initialthermalstate(n,om,r,lambda0,delta,g0,phs,beta)
else    
# initial pure state
psi = DQPT.initialstatequench(n,om,r,lambda0,delta,g0,phs)
psi0 = psi[1]
psi0ad = psi0'
rhoistate = psi0*psi0ad
end
#---------------------------------


#--------building jump operators
adop = diagonalization.anhilation(n)
aop = adop'
Jb = (kp)^(1/2)*aop # bosonic loss
sigmax = (1/2)*diagonalization.sigmax(n)
sigmay = (1/2)*diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
Jq = (kq)^(1/2)*sigmam # qubit relaxation
jpar=[1.0]
jop = [Jb]
# --------------------------------------


function bosonicratefunctiontprime(rho0,Jp,Jq,tmax,n,om,r,lambda,delta,eta,psi,L,tint,acc) 
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 times=(0.0,tmax)
 f(u,p,t) = -im*(HMatrix*u-u*HMatrix) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=0.01,abstol =  acc)
 nt = floor(Int, tmax/tint)
 surprob = []
 tprimelist=[i for i in 0:0.1:2]   
 open("../output/bosonicratefunction_tprime.dat","w") do io
 for tp in tprimelist
   rhotprime = DQPT_thermal.rhot(rhoistate,Jb,Jq,tp,n,om,r,lambda1,delta,g1,phs,L)   
   rhotprime_qo = wigner_eig.buildingrho(rhotprime,n)
   rhotprime_b = ptrace(rhotprime_qo,2)
   println(tp/tprimelist[length(tprimelist)])
   tinst=0.0  
   for i in 1:(nt+1)
     rhot = sol(tinst)  
     rhoqo = wigner_eig.buildingrho(rhot,n)
     rhopt_b = ptrace(rhoqo,2)
     fidb  = tr(rhotprime_b*rhopt_b)
     brf = -(1/r)*log(abs(fidb))  
     println(io,tinst," ",tp," ",brf)
     tinst=tinst+tint
   end
 end
 end
 return "done"
end



#println("...obtaining Lindbladian ")
#HH = diagonalization.hamiltonian(n,om,r,lambda1,delta,g1,phs)
#Lop = diagonalization.Lindbladian(HH,jpar,jop)
#rhot = diagonalization.lindblandDynamics(Lop, rhoistate, tshot)

#wig = wigner_eig.wigner_rhot(rhot,L,r,n)

#println("-> Wigner function obtained via Lindbladian ")    

sp = bosonicratefunctiontprime(rhoistate,Jb,Jq,tmax,n,om,r,lambda1,delta,g1,phs,L,tint,acc)
println(sp)
println("bosonic rate function in file bosonicratefunction_tprime.dat")







end
