module main_dqpts
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/diagonalization.jl")
include("modules/wigner_eig.jl")
include("modules/DQPT.jl")
using .diagonalization
using .DQPT
using .wigner_eig


n=100              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=10.0            # Qubit frequency
lambda0=0.0       # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g0=0.0            # Initial coupling
psi=0.0           # Phase of the Hamiltonian
g1=0.75           # Final coupling
lambda1=0.0       # Final Carrier parameter
nsubint=1000      # Subintervales for integrating the survival probability
nsubint2=400      # Subintervals of real time for estimate the position of zeros
nsubint3=30       # Subintervals of imaginary time for estimate the position of zeros
tmax=10.0         # maximal time for the survival probability
tshot=4.0         # time for the Wigner function
alpha=1.0         # Parameter of the linear combination for the initial state 
ph=0.0            # Phase of the initial state
L=7.5             # Size of the phase space
flag1=0           # (1) for the position of the zeros in the complex plane (0) for skip


name = "output/position_zeros.dat"   # File for saving the position of the zeros

#   Circuits that contains zeros
tcirc=[0.0-0.5*im,0.0+0.5*im,10.0+0.5*im,10.0-0.5*im]
tcircr=[0.0-0.0*im,0.0+0.5*im,10.0+0.5*im,10.0-0.0*im]
tcircl=[0.0-0.5*im,0.0+0.0*im,10.0+0.0*im,10.0-0.5*im]

#------------- Preform calculations---------------------------#



istate = DQPT.initialstatequench(n,om,r,lambda0,delta,g0,psi)

phi0 = alpha^(1/2)*istate[1] + (1-alpha)^(1/2)*exp(im*ph)*istate[2]
hamf = diagonalization.hamiltonian(n,om,r,lambda1,delta,g1,psi)

wigt = DQPT.wigner_rhot(phi0,hamf,L,r,n,tshot)
ctsa = DQPT.amplitud(phi0,tmax,1.0,n,om,r,lambda1,delta,g1,psi,L)


ovl = DQPT.overlapdqpt(phi0,hamf,n)
println("Overlap Done")
rr = DQPT.Nzeros(phi0,tcirc,1.0,n,om,r,lambda1,delta,g1,psi,nsubint)
rrr = DQPT.Nzeros(phi0,tcircr,1.0,n,om,r,lambda1,delta,g1,psi,nsubint)
rrl = DQPT.Nzeros(phi0,tcircl,1.0,n,om,r,lambda1,delta,g1,psi,nsubint)
println("Number of zeros within the full circuit: ",rr)
println("Number of zeros on the right of the real axis: ",rrr)
println("Number of zeros on the left  of the real axis: ",rrl)
println("---------------------------------------------------------")

if flag1==1
 println("-Calculating the position of the zeros in the complex plane-")
 pos = DQPT.PositionsZeros(phi0,tcirc,1.0,n,om,r,lambda1,delta,g1,psi,nsubint2,nsubint3,name)
 println("Number of zeros found in the contour rectangle: ",pos)
 println("Their positions in the complex plane can be found in file position_zeros.dat")
end



end
