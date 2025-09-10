module main_adiabaticramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import adiabatic_ramp



n=100              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=50.0            # Qubit frequency
lambda=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g=(1/2^(1/2))     # Initial coupling
phi=0.0           # Phase of the Hamiltonian
xi=0.05              # modulation amplitude
tau = 4.39        # modulation period
tcycles=250         # number of cycles
L=7.5            # Size of the phase space
k=2               # initial eigenstate
Nf = 1000        # Nfloquet
flagt = 1        # flag for the direction of the perturbation
kf = 1           # target floquet state
acc = 1e-17      # accuracy for the differential equation    



#------------- Perform calculations---------------------------#
psi0 = adiabatic_ramp.initialeigenstateH(n,om,r,lambda,delta,g,phi,k)
println("done")
rhot = adiabatic_ramp.rhoevoladiabatic(psi0,tcycles,n,om,r,lambda,delta,g,phi,xi,tau,Nf,flagt,kf,acc,L)
println("done")
plotwigner = adiabatic_ramp.wignerrhot(rhot,L,r,n)
println("done")

end
