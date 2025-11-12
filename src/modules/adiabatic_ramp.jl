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


function rhot_floquetramp(psi0,tmax,n,om,r,lambda,delta,eta,psi,xi,tau,Nf,flagt,kf,acc,L)
 tfmax = tmax
 tmax = tau*tmax
 psi0t = transpose(conj(psi0))
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 floquet=troterization.troter(n,Nf,r,om,lambda,delta,xi*r,2*pi/tau,eta,psi,flagt)
 vecordered = stat.orderinfvec(n,om,r,lambda,delta,Nf,2*pi/tau,xi*r,eta,psi,flagt)
 evals=eigvals(floquet)
 evs=eigvecs(floquet)
 listvec=[evs[i,vecordered[kf]] for i in 1:length(evals)]
 listvect = transpose(conj(listvec))
 sig = diagonalization.sigmaz(n)
 bc=FockBasis(n)
 ba=SpinBasis(1//2)
 if flagt==2  
     sig = diagonalization.sigmax(n)
 end
 times=(0.0,tmax)
 tint=0.01
 f(u,p,t) = -im*(HMatrix + ((r*xi/2)*((1/2)*(1-cos(pi*t/tmax)))^1)*sig*cos(2*pi*t/tau))*u
 #f(u,p,t) = -im*(HMatrix)*u
 prob = ODEProblem(f,psi0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=tint)
 open("output/fidelities_ramp.dat","w") do io
 open("output/qfi_ramp.dat","w") do io2
 open("output/negativity_ramp.dat","w") do io3
 open("output/magic_ramp.dat","w") do io4        
 #dt=1.0
 dt = tau
 for i in 0:round(Int, tfmax) 
   psit = sol(i*dt)
   #psit = exp(-i*dt*im*HMatrix)*psi0  
   ov1 = psi0t*psit
   ov2 = listvect*psit
   #println("------------>",i*dt," ",abs2(ov1))
   global rhot = psit*transpose(conj(psit))
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   neg  =  wigner_eig.wigner_rhot_neg(rhot,n,L)
   rhotpt = ptrace(rhotqo,2)
   rho_q = ptrace(rhotqo,1)  
   wigner_q = magic.discrete_wigner(rho_q)
   rb = magic.robustness(rho_q)  
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   println(io,i*dt," ",abs2(ov1)," ",abs2(ov2))
   println(io2,i*dt," ",real(qfi)," ",real(qfi2)," ",real(qfi3))
   println(io3,i*dt," ",real(neg))
   println(io4,i*dt," ",real(wigner_q[2]), " ",real(rb))
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
 return rhot
end

function rhot_adiabaticramp(psi0,tmax,n,om,r,lambda0,lambda1,delta,eta,psi,kf,acc,L)
 psi0t = transpose(conj(psi0))
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda1,delta,eta,psi)
 HMatrix0= diagonalization.hamiltonian(n,om,r,lambda0,delta,eta,psi)   
 evals=eigvals(HMatrix)
 evs=eigvecs(HMatrix)
 listvec=[evs[i,kf] for i in 1:length(evals)]
 listvect = transpose(conj(listvec))
 times=(0.0,tmax)
 tint=0.01
 sigmax =  diagonalization.sigmax(n)
 bc=FockBasis(n)
 ba=SpinBasis(1//2)
 #if lambda1>lambda0   
   lam = lambda1 - lambda0
 #else
 #  lam = lambda0 - lambda1
 #end
 f(u,p,t) = -im*(HMatrix0 + (lam/2)*r*(1/2)*(1-cos(pi*t/tmax))*sigmax)*u
 prob = ODEProblem(f,psi0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=tint)
 #dt=1.0
 ddt = 0.1
 open("output/fidelities_ramp.dat","w") do io
 open("output/qfi_ramp.dat","w") do io2
 open("output/negativity_ramp.dat","w") do io3
 open("output/magic_ramp.dat","w") do io4        
 for i in 0:round(Int, tmax/ddt) 
   psit = sol(i*ddt)
   #psit = exp(-i*dt*im*HMatrix)*psi0  
   ov1 = psi0t*psit
   ov2 = listvect*psit
   #println("------------>",i*dt," ",abs2(ov1))
   global rhot = psit*transpose(conj(psit))
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   neg  =  wigner_eig.wigner_rhot_neg(rhot,n,L)
   rhotpt = ptrace(rhotqo,2)
   rho_q = ptrace(rhotqo,1)  
   wigner_q = magic.discrete_wigner(rho_q)
   rb = magic.robustness(rho_q)  
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   println(io,i*ddt," ",abs2(ov1)," ",abs2(ov2))
   println(io2,i*ddt," ",real(qfi)," ",real(qfi2)," ",real(qfi3))
   println(io3,i*ddt," ",real(neg))
   println(io4,i*ddt," ",real(wigner_q[2])," ",real(rb))
 end
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fidelities_ramp.dat to see the fidelities            -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains F(t) from t=0 to ",tmax," in steps of ",ddt," time units             ")
 println("   - Survival probavility                                                                           ")
 println("   - Fidelity between the state and the target state                                                ")   
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_ramp.dat to see the QFI                          -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains QFI(rho(t),A) from t=0 to ",tmax," in steps of ",ddt," time units    ")
 println("   - Phase sensitivity QFI(rho(f),n)                                                                ")
 println("   - Displacement sensitivity QFI(rho(f),p)                                                         ")
 println("   - Displacement sensitivity QFI(rho(f),x)                                                         ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativity_ramp.dat to see the negativity volume       -----------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains delta(t) from t=0 to ",tmax," in steps of ",ddt," time units         ")  
 println("--------------------------------------------------------------------------------------------------- ")
 #psit = sol(tmax)
 #rhot = psit*transpose(conj(psit))
 #rhot = listvec*transpose(conj(listvec))
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
    println("Mana M(rho_q) : ",wigner_q[2])  
 return "done"
end

 
end
