push!(LOAD_PATH, pwd())
module wigner_eig
using LinearAlgebra
using QuantumOptics
using HCubature  # loaded by CalculusWithJulia
using PyPlot
include("stat.jl")
include("diagonalization.jl")
include("troterization.jl")
include("Fisher.jl")
using .stat
using .diagonalization
using .troterization
using .Fisher
#using stat
#using diagonalization
#using troterization
#using Fisher


Qexabs(v) = Qexabs(v...)  # denominator accepts a vector
QWehrl(v) = QWehrl(v...)  # denominator accepts a vector
Psimarginals(v) = Psimarginals(v...)  # denominator accepts a vector


#axis("scaled")


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
  pur=tr(rhopt^2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
  yticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
  #xticks([])
  #yticks([])
  pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
  #colorbar(fraction=0.25,shrink=0.75,ticks=[-0.1,-0.05,0,0.05,0.1])
  #axis("off")
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
  xpm=expect(xop*pop,rhopt)
  pxm=expect(pop*xop,rhopt)
  nexp=expect(adop*aop,rhopt)
  println("central moments of the ",k,"-th stationary state")
  println("x 1 moment (normalized): ",x1m/r^(1/2))
  println("p 1 moment (normalized): ",p1m/r^(1/2))
  println("x 2 moment : ",x2m-x1m^2)
  println("p 2 moment : ",p2m-p1m^2)
  println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  println(" mean photon number n : ",nexp)
  #V11 = x2m-x1m^2
  #V22 = p2m-p1m^2
  #V12 = (1/2)*(xpm+pxm)-x1m*p1m
  #Vmat = [V11 V12; V12 V22]
  #meig = eigvals(Vmat)
  #println("real squeezing parameter: ",(-1/2)*log(2*(real(meig[1]))))
  println("real squeezing parameter: ",real(asinh(nexp^(1/2))))
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  qfi = Fisher.fishern2(rhopt,Nmax)
  qfi2 = Fisher.fisherdisplacementp(rhopt,Nmax)
  qfi3 = Fisher.fisherdisplacementx(rhopt,Nmax)
  println("QFI[rho,n]/(4n) : ", qfi)
  println("QFI[rho,p]/2    : ", qfi2)
  println("QFI[rho,x]/2    : ", qfi3)
  #xy=xycircle((1/(r)^(1/2)),20)
  #scatter(xy[1],xy[2],s=2,color="black")
  #scatter([1],[1],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.001)
  open("output/marginals.dat","w") do io
  global summar3 = 0.0
  for ix in x
    summar1 = 0.0
    summar2 = 0.0
    for ip in x
      summar1 = summar1 + wigner(rhopti, ix, ip)*(2*L/length(x))
      summar2 = summar2 + wigner(rhopti, ip, ix)*(2*L/length(x))
      global summar3 = summar3 + wigner(rhopti, ip, ix)*(2*L/length(x))^2
    end
    println(io,ix," ",round(summar1,digits=8)," ",round(summar2,digits=8))
  end
  end
  #println("********volume:", summar3)
  tight_layout()
  savefig("output/wigner_eigenstate.png")
  println("See output/wigner_eigenstate.png for the Wigner function of the ",k,"-eigenstate of the QRM")
  return [real(den[1]-1.0),pur]
end

function wigner_evolt(Nmax,om,r,lambda,delta,eta,psi,L,phi0,t)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  x = [-L:0.1:L;]
  xm = x/r^(1/2)
  ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
  U = exp(-im*t*ham)
  phit=U*phi0
  listvec=[phit[i] for i in 1:2*Nmax]
  phi = buildingstate(listvec,Nmax)
  rho = dm(phi)
  rhopt = ptrace(rho,2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  xticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
  yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
  #xticks([-2,-1,0,1,2])
  #yticks([-2,-1,0,1,2])
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
  println("central moments of the state at time ",t)
  println("x 1 moment (normalized): ",x1m/r^(1/2))
  println("p 1 moment (normalized): ",p1m/r^(1/2))
  println("x 2 moment : ",x2m-x1m^2)
  println("p 2 moment : ",p2m-p1m^2)
  println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  fotoc = (x2m-x1m^2)^(1/2) + (p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  println("Fotoc : ",real(fotoc))
  xy=xycircle(0.18,20)
  #scatter(xy[1],xy[2],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.0001)
  tight_layout()
  savefig("output/wigner_psit.png")
  println("See output/wigner_psit.png for the Wigner function of the state at time= ", t, " evolving under the time-independent AQRM")
  return [real(den[1]-1.0)]
end

 function wigner_driven(Nmax,om,r,lambda,delta,eta,psi,nu,chi,Nf,k,L,flagt,kl)
   bc=FockBasis(Nmax)
   ba=SpinBasis(1//2)
   x = [-L:0.1:L;]
   xm = x/r^(1/2)
   ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
   floquet=troterization.troter(Nmax,Nf,r,om,lambda,delta,chi,nu,eta,psi,flagt)
   vecordered = stat.orderinfvec(Nmax,om,r,lambda,delta,Nf,nu,chi,eta,psi,flagt)
   evs=eigvecs(floquet)
   listvec=[evs[i,vecordered[k]] for i in 1:2*Nmax]
   evsham = eigvecs(ham)
   listovlp=[]
   listovlpn=[]
   for jk in 1:kl
     klstate=[evsham[i,jk] for i in 1:2*Nmax]
     fockstate=[0.0 for i in 1:2*Nmax]
     fockstate[jk]=1
     ovlp= abs2(transpose(conj(klstate))*listvec)
     ovlpn= abs2(transpose(conj(fockstate))*listvec)
     append!(listovlp,ovlp)
     append!(listovlpn,ovlpn)
   end
   phi = buildingstate(listvec,Nmax)
   rho = dm(phi)
   rhopt = ptrace(rho,2)
   pur=tr(rhopt^2)
   QQ = wigner(rhopt, x, x)
   QQs=transpose(QQ)
   #tick_params(labelsize=30)
   #xticks([-1,0,1])
   #yticks([-1,0,1])
   pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
   global rhopti=rhopt
   den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.001)
   #wentropy=hcubature(QWehrl,(-L,-L),(L,L))
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
   nexp=expect(adop*aop,rhopt)
   println("central moments of the ",k,"-th stationary state")
   println("x 1 moment (normalized): ",x1m/r^(1/2))
   println("p 1 moment (normalized): ",p1m/r^(1/2))
   println("x 2 moment : ",x2m-x1m^2)
   println("p 2 moment : ",p2m-p1m^2)
   println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
   println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
   println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
   println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
   println(" mean photon number n : ",nexp)
   un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
   println("Uncertainity : ",un)
   qfi = Fisher.fishern2(rhopt,Nmax)
   qfi2 = Fisher.fisherdisplacementp(rhopt,Nmax)
   qfi3 = Fisher.fisherdisplacementx(rhopt,Nmax)
   println("QFI[rho,n]/(4n)   : ", qfi)
   println("QFI[rho,p]/2      : ", qfi2)
   println("QFI[rho,x]/2      : ", qfi3)
   tight_layout()
   savefig("output/wigner_eigenstate_f.png")
   println("See output/wigner_eigenstate_f.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
   ham0=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
   eigvecf0=[evs[j,vecordered[k]] for j in 1:2*(Nmax+1)]
   eigvecf0_t = transpose(eigvecf0)
   exp_val = eigvecf0_t*ham0*eigvecf0
   eigvecsH0 = eigvecs(ham0)
   fid = 0.0
   #nk = k
   #for ik in 1:40
   #  statekH0 = [eigvecsH0[j,nk] for j in 1:2*(Nmax+1)]
   #  fidinst = abs2(eigvecf0_t*statekH0)
   #  if fidinst > fid 
   #    fid = fidinst
   #    nk = ik
   #  end
   #end
   return [exp_val,real(den[1]-1.0),pur,listovlp,listovlpn] 
 end

function wigner_evolt_driven(Nmax,om,r,lambda,delta,eta,psi,nu,chi,Nf,L,phi0,pf,flagt)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  x = [-L:0.1:L;]
  xm = x/r^(1/2)
  floquet=troterization.troter(Nmax,Nf,r,om,lambda,delta,chi,nu,eta,psi,flagt)
  phit=phi0
  for i in 1:pf
   phit = floquet*phit
  end
  listvec=[phit[i] for i in 1:2*Nmax]
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
  println("central moments of the state at time ",pf*(2*pi)/nu)
  println("x 1 moment (normalized): ",x1m/r^(1/2))
  println("p 1 moment (normalized): ",p1m/r^(1/2))
  println("x 2 moment : ",x2m-x1m^2)
  println("p 2 moment : ",p2m-p1m^2)
  println("x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  fotoc = (x2m-x1m^2)^(1/2) + (p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  println("Fotoc :",real(fotoc))
  xy=xycircle(0.18,20)
  #scatter(xy[1],xy[2],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L))
  tight_layout()
  savefig("output/wigner_psit_f.png")
  println("See output/wigner_psit_f.png for the Wigner function of the state at time=", pf*(2*pi)/nu," evolving under the Floquet operator of the AQRM")
  return [real(den[1]-1.0)]
end




function wigner_negativities(Nmax,wf,L)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  listvec=[wf[i] for i in 1:2*Nmax]
  phi = buildingstate(listvec,Nmax)
  rho = dm(phi)
  rhopt = ptrace(rho,2)
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.0001)
  return real(den[1]-1.0)
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

function buildingrho(rhoi,Nmax)
rho_ev = eigvals(rhoi)
rho_evecs = eigvecs(rhoi)
bc=FockBasis(Nmax)
ba=SpinBasis(1//2)
psi0= 0.0*fockstate(bc, 0) ⊗ spindown(ba)
rho = dm(psi0)  
for i in 1:length(rho_ev)
   listrho = [rho_evecs[j,i] for j in 1:length(rho_ev)]
   psiinst = buildingstate(listrho,Nmax)
   rhoinst = dm(psiinst)
   #println(rhoinst)
   rho = rho + rho_ev[i]*rhoinst
   #println(rho)
end
return rho
end



function Qexabs(x,y)
   return abs(wigner(rhopti, x, y))
   #return wigner(rhopti, x, y)
   end
   
   

function QWehrl(x,y)
   return qfunc(rhopti, x, y)*log(qfunc(rhopti,x,y))
   #return wigner(rhopti, x, y)
   end

function wigner_rhot(rho,L,r,Nmax)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  x = [-L:0.1:L;]
  xm = x/r^(1/2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  #xy=xycircle((1/(2*r)^(1/2)),20)
  xy=xycircle(0.1,20)
  scatter([0],[0],s=2,color="black")
  pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
  tight_layout()
  savefig("output/wigner_rhot.png")  
  println("See output/wigner_rhot.png for the Wigner function of the state rho(x,p,t)")
  println("The data for figure output/wigner_rhot.png can be found in file output/wigner_rhot_data.dat")
  open("output/wigner_rhot_data.dat","w") do io
  for i in x
    for j in x
      wig = wigner(rhopt, i, j)
      println(io,i/r^(1/2)," ",j/r^(1/2)," ",round(wig,digits=8))
    end
  end
  end
  return "done"
end

function wigner_rhot_neg(rho,Nmax,L)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.0001)
  return real(den[1]-1.0)
end

function wigner_rhot_ZerosCut(rho,Nmax,r)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  maxxp = round(Int,2/0.01)
  int = 0.01
  xplist=[-1+i*int for i in 0:maxxp]
  global count = 0
  #println("-----")
  for i in 1:(length(xplist)-1)
      w1cp = wigner(rhopt, xplist[i] , 0.0)
      w2cp = wigner(rhopt, xplist[i+1], 0.0)
      w1cx = wigner(rhopt, 0.0 , xplist[i])
      w2cx = wigner(rhopt, 0.0, xplist[i+1])
      if (real(w1cp) > 0) && ( real(w2cp) < 0)
         global count = count + 1
	 #println("zero->"," x=",xplist[i]/r^(1/2)," p=",0) 
      end
      if (real(w1cp) < 0) && ( real(w2cp) > 0)
         global count = count + 1
	 #println("zero->"," x=",xplist[i]/r^(1/2)," p=",0) 
      end
      if (real(w1cx) > 0) && ( real(w2cx) < 0)
         global count = count + 1
	 #println("zero->"," p=",xplist[i]/r^(1/2)," x=",0) 
      end
      if (real(w1cx) < 0) && ( real(w2cx) > 0)
         global count = count + 1
	 #println("zero->"," p=",xplist[i]/r^(1/2)," x=",0) 
      end
  end
  return count
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


end

