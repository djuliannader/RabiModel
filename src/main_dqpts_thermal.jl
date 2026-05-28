module main_dqpts_thermal
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/diagonalization.jl")
include("modules/DQPT_thermal.jl")
include("modules/DQPT.jl")
include("modules/wigner_eig.jl")
using .diagonalization
using .DQPT_thermal
using .DQPT
using .wigner_eig



n=60              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=50.0            # Qubit frequency
lambda0=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g0=0              # Initial coxupling
phs=0.0           # Phase of the Hamiltonian
g1=5/4            # Final coupling
lambda1=0.0       # Final Carrier parameter
tshot=7.55         # time for Wigner function
flag1=0           # (1) for the complex time survival probability (0) for skip 
beta=50.0         # temperature
L=10.0             # Size of the phase space
tint = 0.05       # time intervals
acc = 10^(-15)    # accuracy for the diferential equation
kp = 0.1         # photon loss
kq = 0.0          # qubit relaxation    


tmax = 10         # maximal time
tint = 0.05       # time interval
tlist = [i for i in 0:tint:tmax]



#------------- Perform calculations---------------------------#

println("Size of the Fock basis: ",n)
println("Beta:                   ",beta)


#----- Building initial state
# initial thermal state
if beta<100
    rhoistate = diagonalization.initialthermalstate(n,om,r,lambda0,delta,g0,phs,beta)
else    
# initial pure state
psi = DQPT.initialstatequench(n,om,r,lambda0,delta,g0,phs)
psi0 = psi[1]
psi0ad = psi0'
rhoistate = psi0*psi0ad
end
#---------------------------------


#--------building jump operators
aop = diagonalization.anhilation(n)
Jb = (kp)^(1/2)*aop # bosonic loss
sigmax = (1/2)*diagonalization.sigmax(n)
sigmay = (1/2)*diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
Jq = (kq)^(1/2)*sigmam # qubit relaxation
jpar=[1.0]
jop = [Jb]
# --------------------------------------


global wignerft = DQPT_thermal.wignerrhot(rhoistate,Jb,Jq,tshot,n,om,r,lambda1,delta,g1,phs,L)
photondist = DQPT_thermal.photonP(rhoistate,Jb,Jq,tshot,n,om,r,lambda1,delta,g1,phs,L)


#println("...obtaining Lindbladian ")
#HH = diagonalization.hamiltonian(n,om,r,lambda1,delta,g1,phs)
#Lop = diagonalization.Lindbladian(HH,jpar,jop)
#rhot = diagonalization.lindblandDynamics(Lop, rhoistate, tshot)

#wig = wigner_eig.wigner_rhot(rhot,L,r,n)

#println("-> Wigner function obtained via Lindbladian ")    

sp = DQPT_thermal.survivalprobabilityt(rhoistate,Jb,Jq,tlist,n,om,r,lambda1,delta,g1,phs,L,tint,acc)
println(sp)
println("Real-time survival amplitude data in file Loschmidt_amplitud_thermal.dat")





end
