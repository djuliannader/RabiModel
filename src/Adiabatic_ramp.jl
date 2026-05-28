module Adiabatic_ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/adiabatic_ramp.jl")
include("modules/diagonalization.jl")
using .adiabatic_ramp
using .diagonalization



n=50                 # Size of the Fock basis
om=1.0               # Bosonic frequency
r=20.0               # Qubit frequency
lambda0=0.0          # Initial Carrier parameter
lambdaf=0.05         # Final Carrier parameter
delta=0.0            # Parameter (-1,0,1) for (AJC,QRM,JC)
g=1.15               # Coupling strength
phi=0.0              # Phase 
tl = 0.0:2.0:100.0   # list of times for tracking the state
L=7.5                # Size of the phase space
k=1                  # Initial eigenstate of H0
kf = 1               # Target Floquet state
acc = 1e-10          # Accuracy for the differential equation
kb = 0.00            # Bosonic loss
kq = 0.00            # qubit relaxation

#------------- Perform calculations---------------------------#
psi0 = adiabatic_ramp.initialeigenstateH(n,om,r,lambda0,delta,g,phi,k)
println("Building initial state -> Done")
adop = diagonalization.anhilation(n)
aop = adop'
Jb = (kb)^(1/2)*aop # bosonic loss
sigmax = (1/2)*diagonalization.sigmax(n)
sigmay = (1/2)*diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
Jq = (kq)^(1/2)*sigmam # qubit relaxation
rhot = adiabatic_ramp.rhot_adiabaticramp(psi0,tl,n,om,r,lambda0,lambdaf,delta,g,phi,kf,acc,L,Jb,Jq)
println("Solving schrodinger equation -> Done")
plotwigner = adiabatic_ramp.wignerrhot(rhot,L,r,n)
println("Wigner function at the end of the ramp -> Done")
diswigner = adiabatic_ramp.wignerrhot_discrete(rhot,n)



end
