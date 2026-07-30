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




# gate 1
t1 = 4.0
rhot1 = transformations.rhodetcat(Hint1,rho0,t1,n,acc)
# gate 1
t2 = 0.5
rhot2 = transformations.rhodetcat(Hint2,rhot1,t2,n,acc)


t3 = 4.0
g = 1.45
Hqrm  = adop*aop + (r/2)*sigz + g*r^(1/2)*(adop+aop)*(sigx/2) 
rhot3 = transformations.rhodetcat(Hqrm,rhot2,t3,n,acc)

wig = wigner_eig.wigner_rhot(rhot3,L,r,n)


#glist = collect(1.2:0.01:1.8)
#for g in glist
#  Hqrm  = adop*aop + (r/2)*sigz + g*r^(1/2)*(adop+aop)*(sigx/2) 
#  evec = eigvecs(Hqrm)
#  ev =   eigvals(Hqrm)
#  psi = [evec[i,1] for i in 1:length(ev)]
#  psiad = psi'
#  rhoqrm = psi*psiad
#  ov = tr(rhoqrm*rhot2)
#  println("g:",g," ov:", ov)
#end



end
