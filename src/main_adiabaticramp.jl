module main_adiabaticramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/adiabatic_ramp.jl")
using .adiabatic_ramp



n=100              # Size of the Fock basis
om=1.0             # Bosonic frequency
r=50.0             # Qubit frequency
lambda=0.0         # Carrier parameter
delta=0.0          # Parameter (-1,0,1) for (AJC,QRM,JC)
g=(1/2^(1/2))      # Coupling
phi=0.0            # Phase 
xi=0.05            # Modulation amplitude
tau = 8.2          # Modulation period
tcycles=100        # Number of cycles of the ramp
L=10.0             # Size of the phase space
k=1                # Initial eigenstate of H0
Nf = 1000          # Number of subintervals for Trotterization
flagt = 2          # Flag for the time-dependent term (1) for sigma_z (2) for sigma_x
kf = 1             # Target Floquet state
acc = 1e-13        # Accuracy for the differential equation



#------------- Perform calculations---------------------------#
psi0 = adiabatic_ramp.initialeigenstateH(n,om,r,lambda,delta,g,phi,k)
println("Building initial state -> Done")
rhot = adiabatic_ramp.rhoevoladiabatic(psi0,tcycles,n,om,r,lambda,delta,g,phi,xi,tau,Nf,flagt,kf,acc,L)
println("Solving schrodinger equation -> Done")
plotwigner = adiabatic_ramp.wignerrhot(rhot,L,r,n)
println("Wigner function at the end of the ramp -> Done")

end
