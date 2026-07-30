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
g=0.8             # Initial coupling
psi=0.0           # Phase of the Hamiltonia
tint =0.1         # time interval for the survival amplitude
tmax= 10.0        # maximal time for the survival amplitude
tshot=9.55        # time for the Wigner function
alpha=1.0         # Parameter of the linear combination for the initial state 
ph=0.0            # Phase of the initial state
L=10.0             # Size of the phase space
kp = 0.5         # photon loss
kq = 0.0          # qubit relaxation
k = 1             # state sorted according to the real part of the eigenvalue (1) steady state




#------------- Hamiltonian ---------------------------#
Ham = diagonalization.hamiltonian(n,om,r,lambda,delta,g,psi)
#--------building jump operators
aop = diagonalization.anhilation(n)
Jb = (kp)^(1/2)*aop # bosonic loss
sigmax = (1/2)*diagonalization.sigmax(n)
sigmay = (1/2)*diagonalization.sigmay(n)
sigmam = (sigmax-im*sigmay)/2
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
# ---- rho steady -----------
nh = size(Ham, 1)
#global emin = 10
#global imin = 1
#for i in 1:length(LL.values)
#    if abs(real(LL.values[i])) < emin
#        global emin = abs(real(LL.values[i]))
#        global imin = i
#    end
#end
evllist = [abs(real(LL.values[i])) for i in 1:length(LL.values)]
ivl = sortperm(evllist)
rst = LL.vectors[:,ivl[k]]
rhost = reshape(rst, nh, nh)
rhost = 0.5 * (rhost + rhost')
rhost = rhost / real(tr(rhost))
# Wigner function
ww = wigner_eig.wigner_steady(rhost,L,r,n)

#xpos = [0.05*i for i in 0:25]
#for ixpos in xpos
#    cs1 = dynamics.initialcoherent(ixpos*r^(1/2),0.0,0.0,0.0,1.0,n)
#    cs2 = dynamics.initialcoherent(-ixpos*r^(1/2),0.0,0.0,0.0,1.0,n)
#    cs3 = dynamics.initialcoherent(0.0,0.0,0.0,0.0,1.0,n)
#    rhocs1 = cs1*(cs1')
#    rhocs2 = cs2*(cs2')
#    rhocs3 = cs3*(cs3')
#    rhomix = (1/2)*rhocs1 + (1/2)*rhocs2 + 0.0*rhocs3
#    fid = (tr((rhomix^(1/2)*rhost*rhomix^(1/2))^1/2))^2
#    println("xpos: ",ixpos," Fid: ",real(fid))
#end

#xpos = [0.05*i for i in 0:25]
#for ixpos in xpos
#    cs1 = dynamics.initialcoherent(ixpos*r^(1/2),0.0,0.0,0.0,1.0,n)
#    cs2 = dynamics.initialcoherent(-ixpos*r^(1/2),0.0,0.0,0.0,1.0,n)
#    cat = normalize(cs1 + cs2)
#    rhocat = (cat)*(cat')
#    fid = (tr((rhocat^(1/2)*rhost*rhocat^(1/2))^1/2))^2
#    println("xpos: ",ixpos," Fid: ",real(fid))
#end


oxposc = 0.5
cs0 = dynamics.initialcoherent(0.0,0.0,0.0,0.0,1.0,n)
rho0vac = (cs0)*(cs0') 
cs1 = dynamics.initialcoherent(oxposc*r^(1/2),0.0,0.0,0.0,1.0,n)
cs2 = dynamics.initialcoherent(-oxposc*r^(1/2),0.0,0.0,0.0,1.0,n)
cat = normalize(cs1 + cs2)
rho0cat = (cat)*(cat')


#tlist = [2*i for i in 0:500]
tlist = 10 .^ range(-3, 3, length=200)
rhot_mpemba = dynamics.distance_Mpemba(rhost,rho0vac,rho0cat,Ham,Jb,Jq,tlist,1e-16,L,n)

wrhot = wigner_eig.wigner_rhot(rhot_mpemba[2],L,r,n)


end
