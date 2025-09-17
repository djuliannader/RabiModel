module main_dqpts_thermal
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/diagonalization.jl")
include("modules/DQPT_thermal.jl")
using .diagonalization
using .DQPT_thermal



n=80              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=10.0            # Qubit frequency
lambda0=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g0=0              # Initial coupling
psi=0.0           # Phase of the Hamiltonian
g1=3/4            # Final coupling
lambda1=0.0       # Final Carrier parameter
tmax=10.0         # maximal time for the survival probability
tshot=9.5        # time for Wigner function
flag1=0           # (1) for the complex time survival probability (0) for skip 
beta=10.0          # temperature
L=7.5             # Size of the phase space




#------------- Perform calculations---------------------------#

println("Size of the Fock basis: ",n)
println("Beta:                   ",beta)

global rhoistate = DQPT_thermal.initialthermalstate(n,om,r,lambda0,delta,g0,psi,beta)

global wignerft = DQPT_thermal.wignerrhot(rhoistate,tshot,n,om,r,lambda1,delta,g1,psi,L)

sp = DQPT_thermal.survivalprobabilityt(rhoistate,tmax,n,om,r,lambda1,delta,g1,psi,L)
println(sp)
println("Real-time survival amplitude data in file Loschmidt_amplitud_thermal.dat")

if flag1==1 
 sp_c = DQPT_thermal.survivalprobabilityt_ct(rhoistate,[tmax,0.5],n,om,r,lambda1,delta1,eta1,psi1)
 println(sp_c)
end






end
