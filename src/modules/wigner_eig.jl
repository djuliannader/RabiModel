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
include("magic.jl")
using .stat
using .diagonalization
using .troterization
using .Fisher
using .magic
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
  rhor =  listvec*transpose(conj(listvec))    
  rho = dm(phi)
  rhopt = ptrace(rho,2)
  rho_q = ptrace(rho,1)  
  pur=tr(rhopt^2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  xticks([-1,0,1])
  yticks([-1,0,1])
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
  cov=expect(xop*pop+pop*xop,rhopt)  
  println("a) Bosonic mode ")   
  println("Central moments of the ",k,"-th stationary state")
  println("  x 1 moment (normalized): ",x1m/r^(1/2))
  println("  p 1 moment (normalized): ",p1m/r^(1/2))
  println("  x 2 moment : ",x2m-x1m^2)
  println("  p 2 moment : ",p2m-p1m^2)
  println("  x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("  p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("  x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("  p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  println(" Covariance <xp+px>  : ",cov)  
  println("Mean photon number <n> : ",nexp)
  ncuts=10                
  thetal = [i*(pi)/ncuts for i in 0:(ncuts-1)]
  vars = wigner_rhot_variation(rhor,Nmax,r,thetal)
  imvars = argmax(vars)
  #println("Flaaag: ",thetal[imvars]," ",vars[imvars])
  #println("Flaaag: ",thetal)
  #println("Flaaag: ",vars)      
  #V11 = x2m-x1m^2
  #V22 = p2m-p1m^2
  #V12 = (1/2)*(xpm+pxm)-x1m*p1m
  #Vmat = [V11 V12; V12 V22]
  #meig = eigvals(Vmat)
  #println("real squeezing parameter: ",(-1/2)*log(2*(real(meig[1]))))
  println("Real squeezing parameter: ",real(asinh(nexp^(1/2))))
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  qfi = Fisher.fishern2(rhopt,Nmax)
  cfipc = Fisher.cfisherphotoncounting(rhopt,Nmax,[0,pi/2])  
  qfi2 = Fisher.fisherdisplacementp(rhopt,Nmax)
  cfi2 = Fisher.cfisherdisplacement(rhopt,Nmax,0)     
  qfi3 = Fisher.fisherdisplacementx(rhopt,Nmax)
  cfi3 = Fisher.cfisherdisplacement(rhopt,Nmax,pi/2)        
  qfi4 = Fisher.fishersqueezing(rhopt,Nmax)
  println("       Fisher Information (QFI and CFI)  ")
  println("QFI[rho,n]/(4n)            : ", real(qfi))
  println("CFI[rho,n]/(4n)            : ", real(cfipc[1])," photon counting + displacement (beta=",cfipc[2],")")  
  println("QFI[rho,x]/2               : ", real(qfi2))
  println("CFI[rho,x]/2               : ", real(cfi2), " homodyne detection ")
  println("CFI[rho,x]/2               : ", real(cfipc[3])," photon counting + displacement (beta=",cfipc[4],")")           
  println("QFI[rho,p]/2               : ", real(qfi3))
  println("CFI[rho,p]/2               : ", real(cfi3), " homodyne detection ")
  println("CFI[rho,p]/2               : ", real(cfipc[5])," photon counting + displacement (beta=",cfipc[6],")")           
  println("QFI[rho,xp+px]/(4n+2)    : ", real(qfi4))
  #xy=xycircle((1/(r)^(1/2)),20)
  #scatter(xy[1],xy[2],s=2,color="black")
  #scatter([1],[1],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.001)
  minn = 0.0  
  open("output/marginals.dat","w") do io
  global summar3 = 0.0
  for ix in x
    summar1 = 0.0
    summar2 = 0.0
    for ip in x
      summar1 = summar1 + wigner(rhopti, ix, ip)*(2*L/length(x))
      summar2 = summar2 + wigner(rhopti, ip, ix)*(2*L/length(x))
      global summar3 = summar3 + wigner(rhopti, ip, ix)*(2*L/length(x))^2
      winst = wigner(rhopti, ip, ix)
        if winst < minn
            minn = winst
        end
    end
    println(io,ix," ",round(summar1,digits=8)," ",round(summar2,digits=8))
  end
  end
  #println("********volume:", summar3)
  tight_layout()
  cw = wigner_rho(rhopt,L,r,Nmax,"output/wigner_eig_data.dat","output/wignercut_eig_data.dat")
  chusimi = husimi_rho(rhopt,L,r,Nmax,"output/husimi_eig_data.dat")  
  savefig("output/wigner_eigenstate.png")
  println("Bosonic Negativity   delta(rho_b)    : ",real(den[1]-1))
  println("Bosonic Deepest Wigner negativity Wmin: ",minn)  
  println("Purity of the Wigner function        : ",real(pur))
  println( "Bosonic Wehrl entropy               : ",real(chusimi[2]))  
  println("See output/wigner_eigenstate.png for the Wigner function")
     
  mr = magic.robustness(rho_q)
  println("b) Qubit mode*")   
  println("Magic robustness of the qubit mode R(rho_q) : ",mr)
  wigner_q = magic.discrete_wigner(rho_q)
  println("Components of the discrete Wigner function of the qubit mode")
  println("   W00 : ",wigner_q[1][1,1])
  println("   W01 : ",wigner_q[1][1,2])
  println("   W10 : ",wigner_q[1][2,1])
  println("   W11 : ",wigner_q[1][2,2])
  println("Qubit Wigner Negativity delta(rho_q) : ",wigner_q[2])  
  return [real(den[1]-1),pur,minn]
end

function wigner_thermalstate(Nmax,om,r,lambda,delta,eta,psi,L,beta)
  bc=FockBasis(Nmax)
  ba=SpinBasis(1//2)
  x = [-L:0.1:L;]
  xm = x/r^(1/2)
  rho=diagonalization.initialthermalstate(Nmax,om,r,lambda,delta,eta,psi,beta)
  ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
  rhoqo = buildingrho(rho,Nmax)  
  println(" Expectation value of the Hamiltonian: ",real(tr(rho*ham)))
  rhopt = ptrace(rhoqo,2)
  rho_q = ptrace(rhoqo,1)  
  pur=tr(rhopt^2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  xticks([-1,0,1])
  yticks([-1,0,1])
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
  println("a) Bosonic mode ")   
  println("Central moments of the thermal state")
  println("  x 1 moment (normalized): ",x1m/r^(1/2))
  println("  p 1 moment (normalized): ",p1m/r^(1/2))
  println("  x 2 moment : ",x2m-x1m^2)
  println("  p 2 moment : ",p2m-p1m^2)
  println("  x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("  p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("  x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("  p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
  println("Mean photon number n : ",nexp)
  #ncuts=10                
  #thetal = [i*(pi)/ncuts for i in 0:(ncuts-1)]
  #vars = wigner_rhot_variation(rhor,Nmax,r,thetal)
  #imvars = argmax(vars)
  #println("Flaaag: ",thetal[imvars]," ",vars[imvars])
  #println("Flaaag: ",thetal)
  #println("Flaaag: ",vars)      
  #V11 = x2m-x1m^2
  #V22 = p2m-p1m^2
  #V12 = (1/2)*(xpm+pxm)-x1m*p1m
  #Vmat = [V11 V12; V12 V22]
  #meig = eigvals(Vmat)
  #println("real squeezing parameter: ",(-1/2)*log(2*(real(meig[1]))))
  println("Real squeezing parameter: ",real(asinh(nexp^(1/2))))
  un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
  println("Uncertainity : ",un)
  qfi = Fisher.fishern2(rhopt,Nmax)
  cfipc = Fisher.cfisherphotoncounting(rhopt,Nmax,[0,pi/2])  
  qfi2 = Fisher.fisherdisplacementp(rhopt,Nmax)
  cfi2 = Fisher.cfisherdisplacement(rhopt,Nmax,0)     
  qfi3 = Fisher.fisherdisplacementx(rhopt,Nmax)
  cfi3 = Fisher.cfisherdisplacement(rhopt,Nmax,pi/2)        
  qfi4 = Fisher.fishersqueezing(rhopt,Nmax)
  println("       Fisher Information (QFI and CFI)  ")
  println("QFI[rho,n]/(4n)            : ", real(qfi))
  println("CFI[rho,n]/(4n)            : ", real(cfipc[1])," photon counting + displacement (beta=",cfipc[2],")")  
  println("QFI[rho,x]/2               : ", real(qfi2))
  println("CFI[rho,x]/2               : ", real(cfi2), " homodyne detection ")
  println("CFI[rho,x]/2               : ", real(cfipc[3])," photon counting + displacement (beta=",cfipc[4],")")           
  println("QFI[rho,p]/2               : ", real(qfi3))
  println("CFI[rho,p]/2               : ", real(cfi3), " homodyne detection ")
  println("CFI[rho,p]/2               : ", real(cfipc[5])," photon counting + displacement (beta=",cfipc[6],")")           
  println("QFI[rho,xp+px]/(4n+2)    : ", real(qfi4))
  #xy=xycircle((1/(r)^(1/2)),20)
  #scatter(xy[1],xy[2],s=2,color="black")
  #scatter([1],[1],s=2,color="black")
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.001)
  minn = 0.0  
  open("output/marginals.dat","w") do io
  global summar3 = 0.0
  for ix in x
    summar1 = 0.0
    summar2 = 0.0
    for ip in x
      summar1 = summar1 + wigner(rhopti, ix, ip)*(2*L/length(x))
      summar2 = summar2 + wigner(rhopti, ip, ix)*(2*L/length(x))
      global summar3 = summar3 + wigner(rhopti, ip, ix)*(2*L/length(x))^2
      winst = wigner(rhopti, ip, ix)
        if winst < minn
            minn = winst
        end
    end
    println(io,ix," ",round(summar1,digits=8)," ",round(summar2,digits=8))
  end
  end
  #println("********volume:", summar3)
  tight_layout()
  cw = wigner_rho(rhopt,L,r,Nmax,"output/wigner_thermal_data.dat","output/wignercut_themal_data.dat")
  chusimi = husimi_rho(rhopt,L,r,Nmax,"output/husimi_thermal_data.dat")  
  savefig("output/wigner_thermal.png")
  println("Bosonic Negativity   delta(rho_b)    : ",real(den[1]-1))
  println("Bosonic Deepest Wigner negativity Wmin: ",minn)  
  println("Purity of the Wigner function        : ",real(pur))
  println( "Bosonic Wehrl entropy               : ",real(chusimi[2]))  
  println("See output/wigner_thermal.png for the Wigner function")
     
  mr = magic.robustness(rho_q)
  println("b) Qubit mode*")   
  println("Magic robustness of the qubit mode R(rho_q) : ",mr)
  wigner_q = magic.discrete_wigner(rho_q)
  println("Components of the discrete Wigner function of the qubit mode")
  println("   W00 : ",wigner_q[1][1,1])
  println("   W01 : ",wigner_q[1][1,2])
  println("   W10 : ",wigner_q[1][2,1])
  println("   W11 : ",wigner_q[1][2,2])
  println("Qubit Wigner Negativity delta(rho_q) : ",wigner_q[2])  
  return [real(den[1]-1),pur,minn]
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
  println("  x 1 moment (normalized): ",x1m/r^(1/2))
  println("  p 1 moment (normalized): ",p1m/r^(1/2))
  println("  x 2 moment : ",x2m-x1m^2)
  println("  p 2 moment : ",p2m-p1m^2)
  println("  x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
  println("  p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
  println("  x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
  println("  p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
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
  return [real(den[1]-1)]
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
   rho = listvec*(listvec')
   phi = buildingstate(listvec,Nmax)
   rhoqo = dm(phi)
   rhopt = ptrace(rhoqo,2)
   rho_q = ptrace(rhoqo,1)  
   pur=tr(rhopt^2)
   QQ = wigner(rhopt, x, x)
   QQs=transpose(QQ)
   tick_params(labelsize=20)
   xticks([-1,0,1])
   yticks([-1,0,1])
   for jk in 1:kl
     klstate=[evsham[i,jk] for i in 1:2*Nmax]
     fockstatel=[0.0 for i in 1:2*Nmax]
     fockstatel[2*jk-1]=1
     phifock = fockstate(bc,jk-1)
     ovlp= abs2(transpose(conj(klstate))*listvec)
     #ovlpn= abs2(transpose(conj(fockstatel))*listvec)
     ovlpn = expect(rhopt,phifock)
     append!(listovlp,ovlp)
     append!(listovlpn,ovlpn)
   end  
   pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
   global rhopti=rhopt
   #den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.001)
   nneg=wigner_negativities_rhot(Nmax,rho,L)      
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
   cov=expect(xop*pop+pop*xop,rhopt)  
   ham0=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,eta,psi)
   eigvecf0=[evs[j,vecordered[k]] for j in 1:2*(Nmax+1)]
   eigvecf0_t = transpose(eigvecf0)
   exp_val = eigvecf0_t*ham0*eigvecf0
   eigvecsH0 = eigvecs(ham0)
   fid = 0.0
   println("Expectation value      <F_k|H_0|F_k>             : ",real(exp_val))
   println("a) Bosonic mode") 
   println("Central moments of the ",k,"-th stationary state")
   println("  x 1 moment (normalized): ",x1m/r^(1/2))
   println("  p 1 moment (normalized): ",p1m/r^(1/2))
   println("  x 2 moment : ",x2m-x1m^2)
   println("  p 2 moment : ",p2m-p1m^2)
   println("  x 3 moment : ",x3m-3*x2m*x1m+3*x1m^3-x1m^3)
   println("  p 3 moment : ",p3m-3*p2m*p1m+3*p1m^3-p1m^3)
   println("  x 4 moment : ",x4m-4*x3m*x1m+6*x2m*x1m^2-4*x1m^4+x1m^4)
   println("  p 4 moment : ",p4m-4*p3m*p1m+6*p2m*p1m^2-4*p1m^4+p1m^4)
   println(" Covariance <xp+px>  : ",cov)    
   println("Mean photon number n : ",nexp)
   println("Bosonic Negativity   delta(rho_b)    : ",real(2*nneg[2]))  
   println("Purity of the Wigner function        : ",real(pur))
   println(" Root-mean-square Fourier radius of the Wigner function (chi)  : ",real(nneg[3]))  
   un = (x2m-x1m^2)^(1/2)*(p2m-p1m^2)^(1/2)
   println("Uncertainity : ",un)
   qfi = Fisher.fishern2(rhopt,Nmax)
   cfipc = Fisher.cfisherphotoncounting(rhopt,Nmax,[0,pi/2])      
   qfi2 = Fisher.fisherdisplacementp(rhopt,Nmax)
   cfi2 = Fisher.cfisherdisplacement(rhopt,Nmax,0)  
   qfi3 = Fisher.fisherdisplacementx(rhopt,Nmax)
   cfi3 = Fisher.cfisherdisplacement(rhopt,Nmax,pi/2)   
   qfi4 = Fisher.fishersqueezing(rhopt,Nmax)
   println("       Fisher Information (QFI and CFI)  ")
   println("QFI[rho,n]/(4n)            : ", qfi)
   println("CFI[rho,n]/(4n)            : ", real(cfipc[1])," photon counting + displacement (beta=",cfipc[2],")")    
   println("QFI[rho,x]/2               : ", qfi2)
   println("CFI[rho,x]/2               : ", cfi2, " homodyne detection ")
   println("CFI[rho,x]/2            : ", real(cfipc[3])," photon counting (beta=",cfipc[4],")")    
   println("QFI[rho,p]/2               : ", qfi3)
   println("CFI[rho,p]/2               : ", cfi3, " homodyne detection ")
   println("CFI[rho,p]/2               : ", real(cfipc[5])," photon counting + displacement (beta=",cfipc[6],")")             
   println("QFI[rho,xp+px]/(4n+2)    : ", qfi4)
   tight_layout()
   savefig("output/wigner_eigenstate_f.png")
   println("See output/wigner_eigenstate_f.png for the Wigner function of the ",k,"-eigenstate of the Floquet operator")
   cw = wigner_rho(rhopt,L,r,Nmax,"output/wigner_floquet_data.dat","output/wignercut_floquet_data.dat")   
   println("b) Qubit mode*") 
   mr = magic.robustness(rho_q)
   println("Magic robustness of the qubit mode R(rho_q) : ",mr)
   wigner_q = magic.discrete_wigner(rho_q)
   println("Components of the discrete Wigner function of the qubit mode")
   println("   W00 : ",wigner_q[1][1,1])
   println("   W01 : ",wigner_q[1][1,2])
   println("   W10 : ",wigner_q[1][2,1])
   println("   W11 : ",wigner_q[1][2,2])
   println("Qubit Wigner negativity delta(rho_q) : ",wigner_q[2])  
   #nk = k
   #for ik in 1:40
   #  statekH0 = [eigvecsH0[j,nk] for j in 1:2*(Nmax+1)]
   #  fidinst = abs2(eigvecf0_t*statekH0)
   #  if fidinst > fid 
   #    fid = fidinst
   #    nk = ik
   #  end
   #end
   return [exp_val,real(nneg[1]),pur,listovlp,listovlpn] 
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
  println("*Bosonic mode*") 
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
  return [real(den[1]-1)]
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
  return real(den[1]-1)
end

function wigner_negativities_rhot(Nmax,rhot,L)
  rhot_qo = buildingrho(rhot,Nmax)  
  rhot_pt = ptrace(rhot_qo,2)
  xpint = 0.025  
  x = [-L:0.025:L;]
  sumwn = 0.0
  sumwvol = 0.0
  sumchi  = 0.0
  sumdchi  = 0.0  
  for ix in x
    for ip in x
        winst = wigner(rhot_pt, ix, ip)
        winstpx = wigner(rhot_pt, ix+xpint, ip)
        winstpp = wigner(rhot_pt, ix, ip+xpint)
        sumwvol = sumwvol + xpint^2*real(winst)
        dwdx = (winstpx - winst)/xpint 
        dwdp = (winstpp - winst)/xpint    
        if real(winst)<0
            sumwn = sumwn + xpint^2*real(winst)
            sumchi += xpint^2*(dwdx^2 + dwdp^2)
        end
        sumdchi += 4*pi^2*xpint^2*real(winst)^2
    end
  end 
    return [sumwvol,abs(sumwn),real(sumchi/sumdchi)]
end


function wigner_negativities2(Nmax,rhot,rhoin,L,r)
  #bc=FockBasis(Nmax)
  #ba=SpinBasis(1//2)
  #listvec=[wf[i] for i in 1:2*Nmax]
  #phi = buildingstate(listvec,Nmax)
  rhot_qo = buildingrho(rhot,Nmax)  
  global rhot_pt = ptrace(rhot_qo,2)
  rho0_qo = buildingrho(rhoin,Nmax)  
  global rho0_pt = ptrace(rho0_qo,2)  
  #den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.0001)
  x = [-L:0.025:L;]
  minn = 0.0
  gminn = 0.0
  sumwn = 0.0
  sumwvol = 0.0      
  sumwnr = 0.0
  sumip = 0.0
  sumim = 0.0
  sumipl = 0.0
  sumiml = 0.0
  sumipll = 0.0
  sumimll = 0.0  
  ep=10^(-5)  
  for ix in x
    for ip in x
        winst = wigner(rhot_pt, ix, ip)
        w0    = wigner(rho0_pt, ix, ip)
        sumwvol = sumwvol + (0.025)^2*real(winst)
        if real(w0*winst)<(-ep)
            sumim = sumim + (0.025)^2*real(w0*winst)
        else
            sumip = sumip + (0.025)^2*real(w0*winst)
        end 
        if real(winst) < gminn
            gminn = winst
        end
        if real(winst)<(-ep)
                sumwn = sumwn + (0.025)^2*real(winst)
        end
        if (ix^2/r+ip^2/r)<(1/r)
            if real(w0*winst)<(-ep)
              sumiml = sumiml + (0.025)^2*real(w0*winst)
            else
              sumipl = sumipl + (0.025)^2*real(w0*winst)
            end 
            if real(winst)<(-ep)
                sumwnr = sumwnr + (0.025)^2*real(winst)
            end
             if winst < minn
                  minn = winst
             end    
        end    
    end
  end
  for ip in x
     winst = wigner(rhot_pt, 0.0, ip)
     w0    = wigner(rho0_pt, 0.0, ip)
     if real(w0*winst)<(-ep)
            sumimll = sumimll + (0.025)*real(w0*winst)
        else
            sumipll = sumipll + (0.025)*real(w0*winst)
     end  
  end
  return [abs(sumwn),abs(sumwnr),gminn,minn,2*pi*sumip,abs(2*pi*sumim),2*pi*sumipl,abs(2*pi*sumiml),sumwvol,2*pi*sumipll,abs(2*pi*sumimll)]
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
    wabs = abs(wigner(rhopti, x, y))
    ww    = wigner(rhopti, x, y)
   return wabs  
   #return wigner(rhopti, x, y)
   end
   
   

function QWehrl(x,y)
   return qfunc(rhopti, x, y)*log(qfunc(rhopti,x,y))
   #return wigner(rhopti, x, y)
   end

function wigner_rhot(rho,L,r,Nmax)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  x = [-L:0.05:L;]
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

function wigner_steady(rho,L,r,Nmax)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  x = [-L:0.05:L;]
  xm = x/r^(1/2)
  QQ = wigner(rhopt, x, x)
  QQs = transpose(QQ)
  tick_params(labelsize=20)
  #xy=xycircle((1/(2*r)^(1/2)),20)
  xy=xycircle(0.1,20)
  scatter([0],[0],s=2,color="black")
  pcolormesh(xm, xm, QQs, cmap=:bwr,vmin=-0.1,vmax=0.1)
  tight_layout()
  savefig("output/wigner_steady.png")  
  println("See output/wigner_steady.png for the Wigner function of the steady state")
  println("The data for figure output/wigner_steady.png can be found in file output/wigner_steady_data.dat")
  open("output/wigner_steady_data.dat","w") do io
  for i in x
    for j in x
      wig = wigner(rhopt, i, j)
      println(io,i/r^(1/2)," ",j/r^(1/2)," ",round(wig,digits=8))
    end
  end
  end
  return "done"
end


function wigner_rho(rho,L,r,Nmax,name1,name2)
  x = [-L:0.1:L;]
  println("The data for the Wigner function can be found in the file ",name1)
  println("The data for the cuts of the Wigner function can be found in the file ",name2)
  open(name1,"w") do io
  for i in x
    for j in x
      wig = wigner(rho, i, j)
      println(io,i/r^(1/2)," ",j/r^(1/2)," ",round(wig,digits=8))
    end
  end
  end
  open(name2,"w") do io
  for i in x
      wig1 = wigner(rho, i, 0)
      wig2 = wigner(rho, 0, i)
      println(io,i/r^(1/2)," ",round(wig1,digits=8)," ",round(wig2,digits=8))
  end
  end  
  return "done"
end

function husimi_rho(rho,L,r,Nmax,name1)
  x = [-L:0.1:L;]
  hvol  = 0.0 
  wentropy =  0.0 
  println("The data for the Husimi function can be found in the file ",name1)
  open(name1,"w") do io
  for i in x
    for j in x
        hus = qfunc(rho, i, j)
        println(io,i/r^(1/2)," ",j/r^(1/2)," ",round(hus,digits=8))
        hvol = hvol + (1/2)*0.1^2*hus
        wentropy = wentropy + (1/2)*0.1^2*hus*log(hus)  
    end
  end     
  end
  return [wentropy,wentropy] 
end

function wigner_rhot_neg(rho,Nmax,L)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  global rhopti=rhopt
  den=hcubature(Qexabs,(-L,-L),(L,L),rtol=0.0001)
  return real(den[1]-1)
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

function wigner_rhot_ZerosCut2(rho,Nmax,r,thetalist)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  maxxp = round(Int,2/0.01)
  int = 0.01
  inttheta = pi/(2*10)  
  #thetalist  = [i*inttheta for i in 0:10]
  countlist = []  
  for k in 1:length(thetalist)
  theta = thetalist[k]      
  global count = 0
  xlist1 = [(1-int*i)*cos(theta+pi) for i in 0:100]
  xlist2 = [int*i*cos(theta) for i in 1:100]
  xlist=vcat(xlist1,xlist2)
  plist1 = [(1-int*i)*sin(theta+pi) for i in 0:100]
  plist2 = [int*i*sin(theta) for i in 1:100]
  plist=vcat(plist1,plist2)  
  #println("-----")
  for i in 1:(length(xlist)-1)
      w1 = wigner(rhopt, xlist[i] , plist[i])
      w2 = wigner(rhopt, xlist[i+1], plist[i+1])
      if (xlist[i+1]^2/r+plist[i+1]^2/r)^(1/2)<(1/r^(1/2))
      if (real(w1) > 0) && ( real(w2) < 0)
         global count = count + 1
	 #println("zero->"," x=",xplist[i]/r^(1/2)," p=",0) 
      end
      if (real(w1) < 0) && ( real(w2) > 0)
         global count = count + 1
	 #println("zero->"," x=",xplist[i]/r^(1/2)," p=",0) 
      end
      end    
  end
  append!(countlist,count)    
  end    
  return countlist
end

function wigner_rhot_variation(rho,Nmax,r,thetalist)
  rhoqo = buildingrho(rho,Nmax)
  rhopt = ptrace(rhoqo,2)
  maxxp = round(Int,2/0.01)
  int = 0.01
  inttheta = pi/(2*10)  
  #thetalist  = [i*inttheta for i in 0:10]
  countlist = []  
  for k in 1:length(thetalist)
  theta = thetalist[k]      
  global var = 0
  xlist1 = [(1-int*i)*cos(theta+pi) for i in 0:100]
  xlist2 = [int*i*cos(theta) for i in 1:100]
  xlist=vcat(xlist1,xlist2)
  plist1 = [(1-int*i)*sin(theta+pi) for i in 0:100]
  plist2 = [int*i*sin(theta) for i in 1:100]
  plist=vcat(plist1,plist2)  
  #println("-----")
  for i in 1:(length(xlist)-1)
      w1 = wigner(rhopt, xlist[i] , plist[i])
      w2 = wigner(rhopt, xlist[i+1], plist[i+1])
      if (xlist[i+1]^2/r+plist[i+1]^2/r)^(1/2)<(1/r^(1/2))
      #if (real(w1) > 0) && ( real(w2) < 0)
        var = var + abs(w2-w1)
      #end
      #if (real(w1) < 0) && ( real(w2) > 0)
      #  var = var + abs(w2-w1)
      #end
      end    
  end
  append!(countlist,var)    
  end    
  return countlist
end


function wigner_rho_cuts(rho,L,r,Nmax,name1,thetalist)
    rhoqo = buildingrho(rho,Nmax)
    rhopt = ptrace(rhoqo,2)
    maxxp = round(Int,2/0.01)
    int = 0.01
    inttheta = pi/length(thetalist)  
    thetalist  = [i*inttheta for i in 0:length(thetalist)]
    println("The data for the cuts of the Wigner function can be found in the file ",name1)
    open(name1,"w") do io
    for rad in 0:int:L
     wlist=[]
    for i in 1:length(thetalist)
        xinst = (L-rad)*cos(-thetalist[i])
        pinst = (L-rad)*sin(-thetalist[i])        
        wiginst = wigner(rhopt, xinst, pinst)
        append!(wlist,wiginst)    
    end
     println(io,-(L-rad)/r^(1/2)," ",join(wlist," "))
    end
    for rad in int:int:L
    wlist=[]
    for i in 1:length(thetalist)
        xinst = rad*cos(thetalist[i])
        pinst = rad*sin(thetalist[i])        
        wiginst = wigner(rhopt, xinst, pinst)
        append!(wlist,wiginst)    
     end
     println(io,rad/r^(1/2)," ",join(wlist," "))
    end
    end
  return "done"
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

function Nozeroscut(psi,L,N,hbar,xinst,pinst,theta)
       int = 2*L/N
       xlist = [i*sen(theta) for i in -1:int:1]
       plist = [i*cos(theta) for i in -1:int:1]
       count = 0
       for i in 1:(length(xlist)-1)
           w1=wigner.wignerfpoint(wf,L,N,hbar,xlist[i],plist[i])
           w2=wigner.wignerfpoint(wf,L,N,hbar,xlist[i+1],plist[i+1])
       if (real(w1) > 0) && ( real(w2) < 0)
         count = count + 1
       end
       if (real(w1) < 0) && ( real(w2) > 0)
         count = count + 1
       end
       end
       return count
end
