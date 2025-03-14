module main_dqpts_thermal
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
import DQPT_thermal




n=80             # Size of the Fock basis
om=1.0            # Bosonic frequency
r=10.0            # Qubit frequency
lambda0=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g0=0              # Initial coupling
psi=0.0           # Phase of the Hamiltonian
g1=5/4            # Final coupling
lambda1=0.0       # Final Carrier parameter
tmax=10.0         # maximal time for the survival probability
flag1=0           # (1) for the position of the complex time survival probability (0) for skip 
beta=10.0




#------------- Perform calculations---------------------------#

global rhoistate = DQPT_thermal.initialthermalstate(n,om,r,lambda0,delta,g0,psi,beta)

println("Size of the Fock basis: ",n)
println("Beta:                   ",beta)
sp = DQPT_thermal.survivalprobabilityt(rhoistate,tmax,n,om,r,lambda1,delta,g1,psi)
println(sp)
println("Real-time survival amplitud data in file Loschmidt_amplitud_thermal.dat")

if flag1==1 
 sp_c = DQPT_thermal.survivalprobabilityt_ct(rhoistate,[tmax,0.5],n,om,r,lambda1,delta1,eta1,psi1)
 println(sp_c)
end






end
