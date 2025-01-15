push!(LOAD_PATH, pwd())
module wigner_eig
using LinearAlgebra
using QuantumOptics
using HCubature  # loaded by CalculusWithJulia
using PyPlot
import stat
import diagonalization
import troterization


Qexabs(v) = Qexabs(v...)  # denominator accepts a vector

function wigner_eigenstate(Nmax,om,r,lambda,delta,eta,psi,k,L)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  x = [-L:0.1:L;]
  xm = x/r^(1/2)
  ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
  evs=eigvecs(ham)
  listvec=[evs[i,k] for i in 1:2*Nmax]
  phi = buildingstate(listvec,Nmax)
  rho = dm(phi)
  rhopt = ptrace(rho,2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  xticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
  yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
  pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
  #colorbar()
  bc=FockBasis(Nmax)
  adop=create(bc)
  aop = destroy(bc)
  xop=(1/(2)^(1/2))*(aop+adop)
  pop=(im/(2)^(1/2))*(adop-aop)
  x1m=expect(xop,rhopt)
  x2m=expect(xop^2,rhopt)
  x3m=expect(xop^3,rhopt)
  x4m=expect(xop^4,rhopt)
  p1m=expect(pop,rhopt)
  p2m=expect(pop^2,rhopt)
  p3m=expect(pop^3,rhopt)
  p4m=expect(pop^4,rhopt)
  println("central moments of the ",k,"-th stationary state")
  println("x 1 moment (normalized): ",x1m/r^(1/2))
  println("p 1 moment (normalized): ",p1m/r^(1/2))
  println("x 2 moment : ",x2m-x1m^2)
  println("p 2 moment : ",p2m-p1m^2)
  println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  xy=xycircle(0.18,20)
  #scatter(xy[1],xy[2],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L))
  tight_layout()
  savefig("wigner_eigenstate.png")
  println("See wigner_eigenstate.png for the Wigner function of the ",k,"-eigenstate of the QRM")
  return [real(den[1]-1.0)]
end

 function wigner_driven(Nmax,om,r,lambda,delta,eta,psi,nu,chi,Nf,k,L)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   xm = x/r^(1/2)
   floquet=troterization.troter(Nmax,Nf,r,om,lambda,delta,chi,nu,eta,psi)
   vecordered = stat.orderinfvec(Nmax,om,r,lambda,delta,Nf,nu,chi,eta,psi)
   evs=eigvecs(floquet)
   listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
   phi = buildingstate(listvec,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   QQ = wigner(rhopt, x, x)
   QQs=transpose(QQ)
   pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
   global rhopti=rhopt
   den=hcubature(Qexabs,(-L,-L),(L,L))
   bc=FockBasis(Nmax)
   adop=create(bc)
   aop = destroy(bc)
   xop=(1/(2)^(1/2))*(aop+adop)
   pop=(im/(2)^(1/2))*(adop-aop)
   x1m=expect(xop,rhopt)
   x2m=expect(xop^2,rhopt)
   x3m=expect(xop^3,rhopt)
   x4m=expect(xop^4,rhopt)
   p1m=expect(pop,rhopt)
   p2m=expect(pop^2,rhopt)
   p3m=expect(pop^3,rhopt)
   p4m=expect(pop^4,rhopt)
   println("central moments of the ",k,"-th stationary state")
   println("x 1 moment (normalized): ",x1m/r^(1/2))
   println("p 1 moment (normalized): ",p1m/r^(1/2))
   println("x 2 moment : ",x2m-x1m^2)
   println("p 2 moment : ",p2m-p1m^2)
   println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
   println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
   println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
   println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
   un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
   tight_layout()
   savefig("wigner_driven.png")
   println("See wigner_driven.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
   ham0=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
   eigvecf0=[evs[j,vecordered[k]] for j in 1:2*(Nmax+1)]
   eigvecf0_t = transpose(eigvecf0)
   exp_val = eigvecf0_t*ham0*eigvecf0
   return [exp_val,real(den[1]-1.0)] 
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
   tight_layout()
   savefig("wigner_driven2.png")
   println("See wigner_driven2.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
 end

function wigner_driven3(Nmax,om,r,omega,gamma,eta,psi,Nf,k,L,no)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   floquet=troterization.troter3(Nmax,Nf,r,om,gamma,omega,eta,psi,no)
   vecordered = stat.orderinfvec3(Nmax,om,r,gamma,omega,eta,psi,Nf,no)
   #floquet=troterization.troter2(Nmax,Nf,r,om,gamma,omega,eta,psi)
   #vecordered = stat.orderinfvec2(Nmax,om,r,gamma,omega,eta,psi,Nf)
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

function xycircle(r,Npoints)
 xlist=[]
 ylist=[]
 thetaint=2*pi/Npoints
 for i in 1:Npoints
 x=r*cos(thetaint*(i-1))
 y=r*sin(thetaint*(i-1))
 append!(xlist,x) 
 append!(ylist,y)
 end
 return [xlist,ylist]
end

function Qexabs(x,y)
   return abs(wigner(rhopti, x, y))
   end

end

