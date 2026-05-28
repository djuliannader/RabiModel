module main_dqpts
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("modules/diagonalization.jl")
include("modules/wigner_eig.jl")
include("modules/dynamics.jl")
include("modules/transformations.jl")
using .transformations
using .diagonalization
using .dynamics
using .wigner_eig


n=30              # Size of the Fock basis
r=20.0            # Qubit frequency
g=1/2^(1/2)       # coupling
L=10.0             # Size of the phase space
x0 = 0.0   # initial coherent state
p0 = 0.0
theta0 = 0.0
phi0 = 0.0
acc = 10^(-15)



#------------- Perform calculations---------------------------#

psi0 = dynamics.initialcoherent(x0,p0,theta0,phi0,1.0,n)
psi0t = psi0'
rho0 = psi0*psi0t

sigz = diagonalization.sigmaz(n)
sigx = diagonalization.sigmax(n)
sigy = diagonalization.sigmay(n)
aop =  diagonalization.anhilation(n)
adop = aop'
xop = (1/2^(1/2))*(adop+aop)
pop = (im/2^(1/2))*(adop-aop)

Hint1 = pop*sigx
Hint2 = xop*sigy
g=1.2
Hint3 = adop*adop + (r/2)*sigz + (g/2)*r^(1/2)*(adop+aop)*sigx 



# gate 1
t1 = 4.0
rhot1 = transformations.rhodetcat(Hint1,rho0,t1,n,acc)
# gate 1
t2 = 0.5
rhot2 = transformations.rhodetcat(Hint2,rhot1,t2,n,acc)


rhotqo = wigner_eig.wigner_rhot(rhot2,L,r,n)

end
