module Adiabatic_ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/adiabatic_ramp.jl")
using .adiabatic_ramp



n=100              # Size of the Fock basis
om=1.0             # Bosonic frequency
r=10.0             # Qubit frequency
lambda0=0.05       # Initial Carrier parameter
lambdaf=0.0        # Final Carrier parameter
delta=0.0          # Parameter (-1,0,1) for (AJC,QRM,JC)
g=1.2              # Coupling strength
phi=0.0            # Phase 
tramp = 50.0       # Period of the ramp
L=7.5              # Size of the phase space
k=1                # Initial eigenstate of H0
kf = 1             # Target Floquet state
acc = 1e-10        # Accuracy for the differential equation



#------------- Perform calculations---------------------------#
psi0 = adiabatic_ramp.initialeigenstateH(n,om,r,lambda0,delta,g,phi,k)
println("Building initial state -> Done")
rhot = adiabatic_ramp.rhot_adiabaticramp(psi0,tramp,n,om,r,lambda0,lambdaf,delta,g,phi,kf,acc,L)
println("Solving schrodinger equation -> Done")
plotwigner = adiabatic_ramp.wignerrhot(rhot,L,r,n)
println("Wigner function at the end of the ramp -> Done")
diswigner = adiabatic_ramp.wignerrhot_discrete(rhot,n)


end
