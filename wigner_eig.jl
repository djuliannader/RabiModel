push!(LOAD_PATH, pwd())
module wigner_eig
using LinearAlgebra
using QuantumOptics
using PyPlot
import stat
import diagonalization
import troterization


function wigner_eigenstate(Nmax,om,r,lambda,delta,eta,psi,k,L)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  x = [-L:0.1:L;]
  ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
  evs=eigvecs(ham)
  listvec=[evs[i,k] for i in 1:2*Nmax]
  phi = buildingstate(listvec,Nmax)
  rho = dm(phi)
  rhopt = ptrace(rho,2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  pcolormesh(x, x, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
  #title(join("k="*string.(k)))
  tight_layout()
  savefig("wigner_eigenstate.png")
  println("See wigner_eigenstate.png for the Wigner function of the ",k,"-eigenstate of the time-independent Hamiltonian")
end

 function wigner_driven(Nmax,om,r,lambda,delta,eta,psi,nu,chi,Nf,k,L)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   floquet=troterization.troter(Nmax,Nf,r,om,lambda,delta,chi,nu,eta,psi)
   vecordered = stat.orderinfvec(Nmax,om,r,lambda,delta,Nf,nu,chi,eta,psi)
   #println(length(vecordered))
   #println(vecordered)
   evs=eigvecs(floquet)
   listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
   phi = buildingstate(listvec,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   QQ = wigner(rhopt, x, x)
   QQs=transpose(QQ)
   pcolormesh(x, x, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
   #title(join("k="*string.(k)))
   tight_layout()
   savefig("wigner_driven.png")
   println("See wigner_driven.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
 end


function wigner_driven2(Nmax,om,r,omega,gamma,eta,psi,Nf,k,L)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   floquet=troterization.troter2(Nmax,Nf,r,om,gamma,omega,eta,psi)
   vecordered = stat.orderinfvec2(Nmax,om,r,gamma,omega,eta,psi,Nf)
   #println(length(vecordered))
   #println(vecordered)
   evs=eigvecs(floquet)
   listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
   #listvec=[evs[i,1] for i in 1:2*Nmax]
   phi = buildingstate(listvec,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   QQ = wigner(rhopt, x, x)
   QQs = transpose(QQ)
   pcolormesh(x, x, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
   #title(join("k="*string.(k)))
   tight_layout()
   savefig("wigner_driven2.png")
   println("See wigner_driven2.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
 end


 function wigner_drivenqs(Nmax,om,r,lambda,delta,eta,psi,nu,chi,Nf,k,L)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   floquet=troterization.troter(Nmax,Nf,r,om,lambda,delta,chi,nu,eta,psi)
   vecordered = stat.orderinfvecqs(Nmax,om,r,lambda,delta,Nf,nu,chi,eta,psi)
   #println(length(vecordered))
   #println(vecordered)
   evs=eigvecs(floquet)
   listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
   phi = buildingstate(listvec,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   QQ = wigner(rhopt, x, x)
   pcolormesh(x, x, QQ, cmap=:bwr,vmin=-0.1,vmax=0.1)
   #title(join("k="*string.(k)))
   tight_layout()
   savefig("wigner_driven_qs.png")
   println("See wigner_driven.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator sorted according to the quasienergies")
 end
 

function buildingstate(lista,nn)
bc=FockBasis(nn)
ba=SpinBasis(1//2)
psif= 0.0*fockstate(bc, 0) ⊗ spindown(ba)
for i in 0:(nn-1)
    psidown = fockstate(bc, i) ⊗ spindown(ba)
    psiup   = fockstate(bc, i) ⊗ spinup(ba)
    psif = psif + lista[2*i+2]*psiup + lista[2*i+2-1]*psidown 
end
return psif
end


end

