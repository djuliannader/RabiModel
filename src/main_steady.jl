module main_dqpts
push!(LOAD_PATH, pwd())
using LinearAlgebra
using SparseArrays
include("modules/diagonalization.jl")
include("modules/wigner_eig.jl")
include("modules/DQPT.jl")
include("modules/dynamics.jl")
using .diagonalization
using .DQPT
using .wigner_eig
using .dynamics


n=40              # Size of the Fock basis
om=1.0            # Bosonic frequency
r=50.0            # Qubit frequency
lambda=0.0        # Initial Carrier parameter
delta=0.0         # Parameter (-1,0,1) for (AJC,QRM,JC)
g=1.25            # Initial coupling
psi=0.0           # Phase of the Hamiltonia
ph=0.0            # Phase of the initial state
L=10.0            # Size of the phase space
kp = 0.2          # photon loss
kq = 0.0          # qubit relaxation
k = 1             # state sorted according to the real part of the eigenvalue (1) steady state




#------------- Hamiltonian ---------------------------#
Ham = diagonalization.hamiltonian(n,om,r,lambda,delta,g,psi)
aop = diagonalization.anhilation(n)
adop = aop'
sigmaz = diagonalization.sigmaz(n)
sigmax = diagonalization.sigmax(n)
sigmay = diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
#Ham = -1.0*adop*aop + 1*(adop)^2*aop^2 - 3*(aop^2+adop^2)  Kerr
#--------building jump operators
Jb = (kp)^(1/2)*aop # bosonic loss
Jq = (kq)^(1/2)*sigmam # qubit relaxation
jpar=[1.0]
jop = [Jb]
# ----- Linbladian -----------
Linb = diagonalization.Lindbladian(Ham,jpar,jop)
Lop = Matrix(Linb)
LL = eigen(Lop)
open("output/Linbladian_spectrum.dat","w") do io
for i in 1:length(LL.values)
    println(io," ",round(real(LL.values[i]),digits=8)," ",round(imag(LL.values[i])))
end    
end
println("Full Linbladian spectrum in file output/Linbladian_spectrum.dat")
# ---- rho steady -----------
nh = size(Ham, 1)
evllist = [abs(real(LL.values[i])) for i in 1:length(LL.values)]
ivl = sortperm(evllist)
rst = LL.vectors[:,ivl[k]]
rhost = reshape(rst, nh, nh)
rhost = 0.5 * (rhost + rhost')
rhost = rhost / real(tr(rhost))
# Wigner function of the steady state
ww = wigner_eig.wigner_steady(rhost,L,r,n)
trHam= tr(rhost*Ham)
trn= tr(rhost*(adop*aop))
println("tr(rho_s*Ham): ", trHam)
println("tr(rho_s*n): ", trn)



#---  Initial states for Mpemba effect -----#
#cs1 = dynamics.initialcoherent(0.0*r^(1/2),0.15*r^(1/2),0.0,0.0,1.0,n)
#cat = normalize(cs1 + cs2)
#rho0c1 = (cs1)*(cs1')
# time list for relaxation
#tlist = [0,10.0]
#tlist = 10 .^ range(-2, 2, length=100)
#rhot_mpemba = dynamics.distance_Mpemba(rhost,rho0c1,Ham,Jb,Jq,tlist,1e-20,L,n)
# Wigner function at final time
#wrhot = wigner_eig.wigner_rhot(rhot_mpemba,L,r,n)



end
