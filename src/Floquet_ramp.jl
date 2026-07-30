module Floquet_ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/adiabatic_ramp.jl")
include("modules/diagonalization.jl")
using .adiabatic_ramp
using .diagonalization


n=60               # Size of the Fock basis
om=1.0             # Bosonic frequency
r=50.0             # Qubit frequency
lambda=0.0         # Carrier parameter
delta=0.0          # Parameter (-1,0,1) for (AJC,QRM,JC)
g=(1/2^(1/2))      # Coupling
phi=0.0            # Phase 
xi=0.05            # Modulation amplitude
tau = 4.4          # Modulation period
Trcycles=80        # Number of cycles for ramp duration Tr
Tscycles=80        # Number of cycles for the shot of the Wigner function  
L=10.0             # Size of the phase space
k=1                # Initial eigenstate of H0
Nf = 1000          # Number of subintervals for Trotterization
flagt = 1          # Flag for the time-dependent term (1) for sigma_z (2) for sigma_x
kf = 1             # Target Floquet state
acc = 1e-15        # Accuracy for the differential equation
kb  = 0.00         # Bosonic loss
kq  = 0.0          # double photon loss


#

#------------- Perform calculations---------------------------#
psi0 = adiabatic_ramp.initialeigenstateH(n,om,r,lambda,delta,g,phi,k)
println("Building initial state -> Done")
aop = diagonalization.anhilation(n)
Jb = (kb)^(1/2)*aop # bosonic loss
sigmax = (1/2)*diagonalization.sigmax(n)
sigmay = (1/2)*diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
Jq = (kq)^(1/2)*sigmam # qubit relaxation
if kb<10^(-6) && kq < 1e-6
  rhot = adiabatic_ramp.rhot_floquetramp(psi0,Trcycles,Tscycles,n,om,r,lambda,delta,g,phi,xi,tau,Nf,flagt,kf,acc,L)
else
  rhot = adiabatic_ramp.rhot_floquetramp2(psi0,Trcycles,Tscycles,n,om,r,lambda,delta,g,phi,xi,tau,Nf,flagt,kf,acc,L,Jb,Jq)
end
println("Solving schrodinger equation -> Done")
plotwigner = adiabatic_ramp.wignerrhot(rhot[1],L,r,n)
println("Wigner function at the end of the ramp -> Done")
diswigner = adiabatic_ramp.wignerrhot_discrete(rhot[1],n)

#----------  Printing results---------------------------
println("Overlap with the target state at the end of the ramp <F_k|psi(Tramp)>: ",rhot[2])
println("Negativity at the end of the ramp   delta(Tramp): ",2*rhot[3][2])

end
