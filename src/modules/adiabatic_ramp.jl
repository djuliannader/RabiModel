module adiabatic_ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
using QuantumOptics
include("diagonalization.jl")
include("Fisher.jl")
include("troterization.jl")
include("wigner_eig.jl")
include("stat.jl")
include("magic.jl")
using .diagonalization
using .Fisher
using .troterization
using .wigner_eig
using .stat
using .magic
export initialeigenstateH
export rhot_floquetramp
export wignerrhot

function initialeigenstateH(n,om,r,lambda,delta,eta,psi,k)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 evecs=eigvecs(HMatrix)
 evals=eigvals(HMatrix)
 vecn   = [evecs[i,k] for i in 1:length(evals)]
 return vecn
end





function rhot_floquetramp(psi0,Trcycles,Tscycles,n,om,r,lambda,delta,eta,psi,xi,tau,Nf,flagt,kf,acc,L)
 tmax =  tau*Trcycles
 tshot = tau*Tscycles
 tlist = [tau*i for i in 0:Trcycles]
 times = (0.0,tmax)      
 psi0t = transpose(conj(psi0))
 rho0 = psi0*(psi0')
 rho0_qo = stat.buildingrho(rho0,n)
 rho0_pt = ptrace(rho0_qo,2)   
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 floquet=troterization.troter(n,Nf,r,om,lambda,delta,xi*r,2*pi/tau,eta,psi,flagt)
 vecordered = stat.orderinfvec(n,om,r,lambda,delta,Nf,2*pi/tau,xi*r,eta,psi,flagt)
 evals=eigvals(floquet)
 evs=eigvecs(floquet)
 listvec=[evs[i,vecordered[kf]] for i in 1:length(evals)]
 listvect = transpose(conj(listvec))
 rhofl = listvec*(listvect)
 rhofl_qo = stat.buildingrho(rhofl,n)
 rhofl_pt = ptrace(rhofl_qo,2)      
 sig = diagonalization.sigmaz(n)
 bc=FockBasis(n)
 ba=SpinBasis(1//2)
 adop=create(bc)
 aop = destroy(bc)
 xop=(1/(2)^(1/2))*(aop+adop)
 pop=(im/(2)^(1/2))*(adop-aop)   
 ov2=0
 neg=0
 man=0.0
 rob = 0.0
 qfi = 0.0
 qfi2 = 0.0
 qfi3 = 0.0   
 if flagt==2  
     sig = diagonalization.sigmax(n)
 end
 tint=0.1
 f(u,p,t) = -im*(HMatrix + ((r*xi/2)*((1/2)*(1-cos(pi*t/tmax)))^1)*sig*cos(2*pi*t/tau))*u
 prob = ODEProblem(f,psi0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff], dt=0.1, saveat = tlist, save_everystep = false)
 open("output/fidelities_ramp.dat","w") do io
 open("output/qfi_ramp.dat","w") do io2
 open("output/negativity_ramp.dat","w") do io3
 open("output/magic_ramp.dat","w") do io4        
 #dt=1.0
 dt = tau
 for i in eachindex(sol.t)
   println(sol.t[i])              
   psit = sol.u[i]
   psitt = psit'
   rhot = psit*psitt
   rhot_qo = stat.buildingrho(rhot,n)
   rhotpt =  ptrace(rhot_qo,2)    
   nfac = psitt*psit  
   psit = psit/nfac^(1/2)
   cor = stat.negativitycorr(rhofl,rhot,n,L,r)    
   #println("------------>",i*dt," ",abs2(ov1))
   rhot = psit*transpose(conj(psit))
   if abs(sol.t[i]-tshot)<10^(-6)
     global rhots = rhot
   end    
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   rhotpt = ptrace(rhotqo,2)
   #ov1 = fidelity(rhotpt,rho0_pt)
   #ov2 = fidelity(rhotpt,rhofl_pt)
   ov1 = psi0t*rhot*psi0
   ov2 = listvect*rhot*listvec  
   neg  =  wigner_eig.wigner_negativities_rhot(n,rhot,L)  
   rho_q = ptrace(rhotqo,1)  
   wigner_q = magic.discrete_wigner(rho_q)
   rb = magic.robustness(rho_q)  
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   man = real(wigner_q[2])
   rob = real(rb)
   #-----------Reference Gaussian state
   x1m = real(expect(xop, rhotpt))
   x2m = real(expect(xop^2, rhotpt))
   p1m = real(expect(pop, rhotpt))
   p2m = real(expect(pop^2, rhotpt))
   varx = x2m - x1m^2
   varp = p2m - p1m^2
   xp_sym = real(expect(xop*pop + pop*xop, rhotpt))
   C = xp_sym/2 - x1m*p1m
   detV = varx*varp - C^2
   nu = sqrt(detV)
   nt = nu -0.5
   ps = 0.5 * atan(2*real(C), real(varx) - real(varp))
   rxi = 0.25 * log(real(varx)/real(varp))
   xi = rxi*exp(im*ps)  
   alpha = (1/(2)^(1/2))*(x1m+im*p1m)  
   Sqop = exp((1/2)*xi*adop^2 - (1/2)*conj(xi)*aop^2)
   Dop = exp(alpha*adop - conj(alpha)*aop)
   #  println("----Flag----")
   #  println("nu:",nu)
   #  println("xi:",xi)
   #  println("xav:",real(x1m)/r^(1/2))
   #  println("pav:",real(p1m)/r^(1/2))
   state0 = fockstate(bc,0)
   rho0 = dm(state0)
   rho_gr = rho0
   rho_gr = (nt^0/(1+nt))*rho0
   for nk in 1:n
       staten = fockstate(bc,nk)
       rhon = dm(staten)
       rho_gr = rho_gr + (nt^(nk)/(1+nt)^(nk+1))*rhon
   end
   rho_gr = Sqop*rho_gr*(Sqop')
   rho_gr = Dop*rho_gr*(Dop')
   qfi_gr = Fisher.fishern2(rho_gr,n)
   qfi2_gr = Fisher.fisherdisplacementp(rho_gr,n)
   qfi3_gr = Fisher.fisherdisplacementx(rho_gr,n)
   # println("qfi:",qfi_gr) 
   #println("------------")    
   #-----------  
   println(io,sol.t[i]," ",real(ov1)," ",real(ov2)," ",real(cor[1]))
   println(io2,sol.t[i]," ",real(qfi)," ",real(qfi2)," ",real(qfi3)," ",real(qfi_gr)," ",real(qfi2_gr)," ",real(qfi3_gr))
   println(io3,sol.t[i]," ",(join(neg, " ")))
   println(io4,sol.t[i]," ",man, " ",rob)
 end
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fidelities_ramp.dat to see the fidelities            -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains F(t) from t=0 to ",tmax," in steps of ",tau," time units             ")
 println("   - Survival probavility                                                                           ")
 println("   - Fidelity between the state and the target state                                                ")   
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_ramp.dat to see the QFI                          -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains QFI(rho(t),A) from t=0 to ",tmax," in steps of ",tau," time units    ")
 println("   - Phase sensitivity QFI(rho(f),n)                                                                ")
 println("   - Displacement sensitivity QFI(rho(f),p)                                                         ")
 println("   - Displacement sensitivity QFI(rho(f),x)                                                         ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativity_ramp.dat to see the negativity volume       -----------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains delta(t) from t=0 to ",tmax," in steps of ",tau," time units         ")  
 println("--------------------------------------------------------------------------------------------------- ")
 println("Final fidelity to target state : ", abs2(ov2))
 println("Final QFI phase sensing (n)    : ", real(qfi))
 println("Final QFI displacement (x)     : ", real(qfi2))
 println("Final QFI displacement (p)     : ", real(qfi3))
 println("Final bosonic negativity       : ", neg)   
 println("Final qubit negativity         : ", man)
 println("Final qubit magic              : ", rob)
 return [rhots,abs2(ov2),real(neg)]
end


function rhot_floquetramp2(psi0,Tcycles,Tsc,n,om,r,lambda,delta,eta,psi,xi,tau,Nf,flagt,kf,acc,L,Jp,Jq)
 tmax =  tau*Tcycles
 tshot = tau*Tsc  
 tlist = [tau*i for i in 0:Tcycles]
 times = (0.0,tmax)      
 psi0t = transpose(conj(psi0))
 rho0 = psi0*(psi0')
 rho0_qo = stat.buildingrho(rho0,n)
 rho0_pt = ptrace(rho0_qo,2)   
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 floquet=troterization.troter(n,Nf,r,om,lambda,delta,xi*r,2*pi/tau,eta,psi,flagt)
 vecordered = stat.orderinfvec(n,om,r,lambda,delta,Nf,2*pi/tau,xi*r,eta,psi,flagt)
 evals=eigvals(floquet)
 evs=eigvecs(floquet)
 listvec=[evs[i,vecordered[kf]] for i in 1:length(evals)]
 listvect = transpose(conj(listvec))
 rhofl = listvec*(listvect)
 rhofl_qo = stat.buildingrho(rhofl,n)
 rhofl_pt = ptrace(rhofl_qo,2)      
 sig = diagonalization.sigmaz(n)
 bc=FockBasis(n)
 ba=SpinBasis(1//2)
 adop=create(bc)
 aop = destroy(bc)
 xop=(1/(2)^(1/2))*(aop+adop)
 pop=(im/(2)^(1/2))*(adop-aop)   
 ov2 = 0.0 
 bnegf = 0.0 
 qrobf = 0.0        
 qnegf = 0.0
 qfi = 0.0
 qfi2 = 0.0
 qfi3 = 0.0
 if flagt==2  
     sig = diagonalization.sigmax(n)
 end
 tint=0.1
    f(u,p,t) = -im*(HMatrix + ((r*xi/2)*((1/2)*(1-cos(pi*t/tmax)))^1)*sig*cos(2*pi*t/tau))*u + im*u*(HMatrix + ((r*xi/2)*((1/2)*(1-cos(pi*t/tmax)))^1)*sig*cos(2*pi*t/tau)) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff], dt=0.1, saveat = tlist, save_everystep = false)
 println("-> PDE solved ")      
 open("output/fidelities_ramp.dat","w") do io
 open("output/qfi_ramp.dat","w") do io2
 open("output/negativity_ramp.dat","w") do io3
 open("output/magic_ramp.dat","w") do io4        
 #dt=1.0
 dt = tau
 for i in eachindex(sol.t)
   println(sol.t[i])              
   global rhot = sol.u[i]
   rhot = (rhot + rhot')/2
   fac = tr(rhot)
   rhot = rhot/fac
   if abs(sol.t[i]-tshot)<10^(-8)
       global rhots = rhot
   end       
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   rhotpt = ptrace(rhotqo,2)  
   #ov1 = fidelity(rhotpt,rho0_pt)
   #ov2 = fidelity(rhotpt,rhofl_pt)
   ov1 = psi0t*rhot*psi0
   ov2 = listvect*rhot*listvec
   cor = stat.negativitycorr(rhofl,rhot,n,L,r)    
   global neg  =  wigner_eig.wigner_negativities_rhot(n,rhot,L)
   rho_q = ptrace(rhotqo,1)
   qnegf = magic.discrete_wigner(rho_q)
   qrobf = magic.robustness(rho_q)
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   #-----------Reference Gaussian state
   x1m = real(expect(xop, rhotpt))
   x2m = real(expect(xop^2, rhotpt))
   p1m = real(expect(pop, rhotpt))
   p2m = real(expect(pop^2, rhotpt))
   varx = x2m - x1m^2
   varp = p2m - p1m^2
   xp_sym = real(expect(xop*pop + pop*xop, rhotpt))
   C = xp_sym/2 - x1m*p1m
   detV = varx*varp - C^2
   nu = sqrt(detV)
   nt = nu -0.5
   ps = 0.5 * atan(2*real(C), real(varx) - real(varp))
   rxi = 0.25 * log(real(varx)/real(varp))
   xi = rxi*exp(im*ps)  
   alpha = (1/(2)^(1/2))*(x1m+im*p1m)  
   Sqop = exp((1/2)*xi*adop^2 - (1/2)*conj(xi)*aop^2)
   Dop = exp(alpha*adop - conj(alpha)*aop)
   state0 = fockstate(bc,0)
   rho0 = dm(state0)
   rho_gr = rho0
   rho_gr = (nt^0/(1+nt))*rho0
   for nk in 1:n
       staten = fockstate(bc,nk)
       rhon = dm(staten)
       rho_gr = rho_gr + (nt^(nk)/(1+nt)^(nk+1))*rhon
   end
   rho_gr = Sqop*rho_gr*(Sqop')
   rho_gr = Dop*rho_gr*(Dop')
   qfi_gr = Fisher.fishern2(rho_gr,n)
   qfi2_gr = Fisher.fisherdisplacementp(rho_gr,n)
   qfi3_gr = Fisher.fisherdisplacementx(rho_gr,n)
   #-----------    
   println(io,sol.t[i]," ",abs(ov1)," ",abs(ov2)," ",real(cor[1]))
   println(io2,sol.t[i]," ",real(qfi)," ",real(qfi2)," ",real(qfi3)," ",real(qfi_gr)," ",real(qfi2_gr)," ",real(qfi3_gr))
   println(io3,sol.t[i]," ",(join(neg, " ")))
   println(io4,sol.t[i]," ",real(qnegf[2]), " ",real(qrobf))
 end
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fidelities_ramp.dat to see the fidelities            -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains F(t) from t=0 to ",tmax," in steps of ",tau," time units             ")
 println("   - Survival probavility                                                                           ")
 println("   - Fidelity between the state and the target state                                                ")   
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_ramp.dat to see the QFI                          -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains QFI(rho(t),A) from t=0 to ",tmax," in steps of ",tau," time units    ")
 println("   - Phase sensitivity QFI(rho(f),n)                                                                ")
 println("   - Displacement sensitivity QFI(rho(f),p)                                                         ")
 println("   - Displacement sensitivity QFI(rho(f),x)                                                         ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativity_ramp.dat to see the negativity volume       -----------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains delta(t) from t=0 to ",tmax," in steps of ",tau," time units         ")  
 println("--------------------------------------------------------------------------------------------------- ")
 #psit = sol(tmax)
 #rhot = psit*transpose(conj(psit))
 #rhot = listvec*transpose(conj(listvec))
 println("Final fidelity to target state : ", real(ov2))
 println("Final QFI phase sensing (n)    : ", real(qfi))
 println("Final QFI displacement (x)     : ", real(qfi2))
 println("Final QFI displacement (p)     : ", real(qfi3))   
 println("Final bosonic negativity       : ", neg[1])
 println("Final qubit negativity         : ", real(qnegf[2]))
 println("Final qubit magic              : ", real(qrobf))
 return [rhots,abs2(ov2),neg]
end


function rhot_adiabaticramp(HMatrix0,rho0,rhotar,tlist,n,om,r,lambda0,lambda1,g0,g1,acc,L,Jp,Jq)   
 tmax = tlist[length(tlist)]
 times = (tlist[1],tmax)   
 sigmax =  diagonalization.sigmax(n)
 aop = diagonalization.anhilation(n)
 adop = aop'
 bc=FockBasis(n)
 ba=SpinBasis(1//2)
 adopqo=create(bc)
 aopqo = destroy(bc)
 xop=(1/(2)^(1/2))*(aopqo+adopqo)
 pop=(im/(2)^(1/2))*(adopqo-aopqo)   
 lam = lambda1 - lambda0
 gdif = g1 -g0
    f(u,p,t) = -im*(HMatrix0 + (sigmax*(lam*r/2) + gdif*r^(1/2)*(aop+adop)*sigmax/2)*(1/2)*(1-cos(pi*t/tmax)))*u + im*u*(HMatrix0 + (sigmax*(lam*r/2) + gdif*r^(1/2)*(aop+adop)*sigmax/2)*(1/2)*(1-cos(pi*t/tmax))) + (Jp*u*(Jp') - ((Jp')*Jp*u + u*(Jp')*Jp)/2) + (Jq*u*(Jq') - ((Jq')*Jq*u + u*(Jq')*Jq)/2)
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=0.1,saveat = tlist,
             save_everystep = false)
 ov2 = 0.0 
 bnegf = 0.0 
 robf = 0.0        
 qnegf = 0.0
 wminf = 0.0
 impf = 0.0
 varpf = 0.0  
 println("-> PDE solved ")   
 open("output/fidelities_ramp.dat","w") do io
 open("output/qfi_ramp.dat","w") do io2       
 open("output/negativity_ramp.dat","w") do io3
 open("output/magic_ramp.dat","w") do io4
 for i in eachindex(sol.t)
   println(sol.t[i])        
   global rhot = sol.u[i]
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   global rhotpt = ptrace(rhotqo,2)  
   ov1 = real(tr(sqrt(sqrt(Hermitian(rho0)) * rhot * sqrt(Hermitian(rho0))))^2)
   ov2 = real(tr(sqrt(sqrt(Hermitian(rhotar)) * rhot * sqrt(Hermitian(rhotar))))^2)
   cor = stat.negativitycorr(rhotar,rhot,n,L,r)      
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   global neg  =  wigner_eig.wigner_negativities_rhot(n,rhot,L)
   bnegf = 2*neg[2]
   wminf = neg[4]
   impf = tr(rhotpt^2)  
   rho_q = ptrace(rhotqo,1)
   wigner_q = magic.discrete_wigner(rho_q)
   qnegf = real(wigner_q[2])  
   rb = magic.robustness(rho_q)
   robf = real(rb)
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   #-----------Reference Gaussian state
   x1m = real(expect(xop, rhotpt))
   x2m = real(expect(xop^2, rhotpt))
   p1m = real(expect(pop, rhotpt))
   p2m = real(expect(pop^2, rhotpt))
   varx = x2m - x1m^2
   varp = p2m - p1m^2
   varpf = varp    
   xp_sym = real(expect(xop*pop + pop*xop, rhotpt))
   C = xp_sym/2 - x1m*p1m
   detV = varx*varp - C^2
   nu = sqrt(detV)
   nt = nu -0.5
   ps = 0.5 * atan(2*real(C), real(varx) - real(varp))
   rxi = 0.25 * log(real(varx)/real(varp))
   xi = rxi*exp(im*ps)
   alpha = (1/(2)^(1/2))*(x1m+im*p1m)  
   Sqop = exp((1/2)*xi*adopqo^2 - (1/2)*conj(xi)*aopqo^2)
   Dop = exp(alpha*adopqo - conj(alpha)*aopqo)
   state0 = fockstate(bc,0)
   rho00 = dm(state0)
   rho_gr = rho00
   rho_gr = (nt^0/(1+nt))*rho00
   for nk in 1:n
       staten = fockstate(bc,nk)
       rhon = dm(staten)
       rho_gr = rho_gr + (nt^(nk)/(1+nt)^(nk+1))*rhon
   end
   rho_gr = Sqop*rho_gr*(Sqop')
   rho_gr = Dop*rho_gr*(Dop')
   qfi_gr = Fisher.fishern2(rho_gr,n)
   qfi2_gr = Fisher.fisherdisplacementp(rho_gr,n)
   qfi3_gr = Fisher.fisherdisplacementx(rho_gr,n)  
   #-----------   
   println(io,sol.t[i]," ",real(ov1)," ",real(ov2)," ",real(cor[1]))
   println(io2,sol.t[i]," ",real(qfi)," ",real(qfi2)," ",real(qfi3)," ",real(qfi_gr)," ",real(qfi2_gr)," ",real(qfi3_gr))  
   println(io3,sol.t[i]," ",(join(neg, " ")))
   println(io4,sol.t[i]," ",qnegf," ",robf)
 end
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fidelities_ramp.dat to see the fidelities            -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains F(t) from t=0 to ",tmax," in steps of ",tlist[2]-tlist[1]," time units             ")
 println("   - Survival probavility                                                                           ")
 println("   - Fidelity between the state and the target state                                                ")   
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_ramp.dat to see the QFI                          -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains QFI(rho(t),A) from t=0 to ",tmax," in steps of ",tlist[2]-tlist[1]," time units    ")
 println("   - Phase sensitivity QFI(rho(f),n)                                                                ")
 println("   - Displacement sensitivity QFI(rho(f),p)                                                         ")
 println("   - Displacement sensitivity QFI(rho(f),x)                                                         ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativity_ramp.dat to see the negativity volume       -----------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains delta(t) from t=0 to ",tmax," in steps of ",tlist[2]-tlist[1]," time units         ")  
 println("--------------------------------------------------------------------------------------------------- ")
 println("Final fidelity to target state : ", real(ov2))
 println("Final bosonic negativity       : ", bnegf)
 println("deepest value of Wigner function  : ", wminf)
 println("Oscillator impurity               : ", real(impf))
 println("Final Var(P)               : ", varpf)   
 cf = Fisher.cfisherphotoncounting(rhotpt,n,[0,pi/2])
 println("CFI                               : ", real(cf[3])," photon counting + displacement (beta=",cf[4],")")            
 println("Final qubit negativity         : ", qnegf)
 println("Final qubit magic              : ", robf)
 open("output/marginals.dat","w") do io
 x = [-L:0.1:L;]       
  for ix in x
    summar1 = 0.0
    summar2 = 0.0
    for ip in x
      summar1 = summar1 + wigner(rhotpt, ix, ip)*(2*L/length(x))
      summar2 = summar2 + wigner(rhotpt, ip, ix)*(2*L/length(x))
    end
    println(io,ix," ",round(summar1,digits=8)," ",round(summar2,digits=8))
  end
  end
 return rhot
end



function wignerrhot(rhot,L,r,n)
 wig = wigner_eig.wigner_rhot(rhot,L,r,n)
 return "done"
end

function wignerrhot_discrete(rhot,n)
    bc=FockBasis(n)
    ba=SpinBasis(1//2)
    rhotqo =  wigner_eig.buildingrho(rhot,n)
    rho_q = ptrace(rhotqo,1)  
    wigner_q = magic.discrete_wigner(rho_q)
    println("Components of the discrete Wigner function of the qubit mode")
    println("   W00 : ",wigner_q[1][1,1])
    println("   W01 : ",wigner_q[1][1,2])
    println("   W10 : ",wigner_q[1][2,1])
    println("   W11 : ",wigner_q[1][2,2])
    println("Qubit Wigner-negativity delta(rho_q) : ", wigner_q[2])  
 return "done"
end

 
end
