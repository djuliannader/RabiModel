module Adiabatic_ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/adiabatic_ramp.jl")
include("modules/diagonalization.jl")
using .adiabatic_ramp
using .diagonalization



n=60                 # Size of the Fock basis
om=1.0               # Bosonic frequency
r=50.0               # Qubit frequency
delta=0.0            # Parameter (-1,0,1) for (AJC,QRM,JC)
phi=0.0              # Phase 
lambda0=0.0          # Initial assymetry bias
lambda1=0.0          # Final assymetry bias
g0=0.0               # Initial Coupling strength
g1=1.135               # Final Coupling strength
tl = 0.0:5.0:200.0   # list of times for tracking the state
L=10.0               # Size of the phase space
k=1                  # Initial eigenstate of H0
kf = 1               # Target Floquet state
acc = 1e-22          # Accuracy for the differential equation
kb = 0.00100             # Bosonic loss
kq = 0.00            # qubit relaxation

#------------- Perform calculations---------------------------#
aop = diagonalization.anhilation(n)
adop = aop'
Jb = (kb)^(1/2)*aop # bosonic loss
sigmaz = diagonalization.sigmaz(n)
sigmax = diagonalization.sigmax(n)
sigmay = diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
Jq = (kq)^(1/2)*sigmam # qubit relaxation
H0 = adop*aop + (r/2)*sigmaz + g0*r^(1/2)*(adop + aop)*sigmax/2  + lambda0*r*sigmax/2
H1 = adop*aop + (r/2)*sigmaz + g1*r^(1/2)*(adop + aop)*sigmax/2  + lambda1*r*sigmax/2
evecs0 = eigvecs(H0)
psi0=evecs0[:,1]
rho0=psi0*(psi0')
evecs1 = eigvecs(H1)
psi1=evecs1[:,1]
rho1=psi1*(psi1')
println("Obtaining initial state -> Done")
rhot = adiabatic_ramp.rhot_adiabaticramp(H0,rho0,rho1,tl,n,om,r,lambda0,lambda1,g0,g1,acc,L,Jb,Jq)
println("Solving the master equation -> Done")
plotwigner = adiabatic_ramp.wignerrhot(rhot,L,r,n)
println("Wigner function at the end of the ramp -> Done")
diswigner = adiabatic_ramp.wignerrhot_discrete(rhot,n)



end
