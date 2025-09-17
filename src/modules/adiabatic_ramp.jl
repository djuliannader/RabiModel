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
using .diagonalization
using .Fisher
using .troterization
using .wigner_eig
using .stat
export initialeigenstateH
export rhoevoladiabatic
export wignerrhot

function initialeigenstateH(n,om,r,lambda,delta,eta,psi,k)
 HMatrix= diagonalization.hamiltonian(n,om,r,lambda,delta,eta,psi)
 evecs=eigvecs(HMatrix)
 evals=eigvals(HMatrix)
 vecn   = [evecs[i,k] for i in 1:length(evals)]
 return vecn
end


function rhoevoladiabatic(psi0,tmax,n,om,r,lambda,delta,eta,psi,xi,tau,Nf,flagt,kf,acc,L)
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
 if flagt==2  
     sig = diagonalization.sigmax(n)
 end
 times=(0.0,tmax)
 tint=0.01
 f(u,p,t) = -im*(HMatrix + ((r*xi/2)*((1/2)*(1-cos(pi*t/tmax)))^1)*sig*cos(2*pi*t/tau))*u
 #f(u,p,t) = -im*(HMatrix + (r*xi/(2))*sigz*cos(2*pi*t/tau))*u
 prob = ODEProblem(f,psi0,times)
 #prob = ODEProblem(f,listvec,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=tint)
 open("output/fidelities_adiabatic.dat","w") do io
 open("output/qfi_adiabatic.dat","w") do io2
 open("output/negativity_adiabatic.dat","w") do io3
 #dt=1.0
 dt = tau
 for i in 0:round(Int, tfmax) 
   psit = sol(i*dt)
   ov1 = psi0t*psit
   ov2 = listvect*psit
   global rhot = psit*transpose(conj(psit))
   rhotqo =  wigner_eig.buildingrho(rhot,n)
   neg  =  wigner_eig.wigner_rhot_neg(rhot,n,L)
   rhotpt = ptrace(rhotqo,2)
   qfi = Fisher.fishern2(rhotpt,n)
   qfi2 = Fisher.fisherdisplacementp(rhotpt,n)
   qfi3 = Fisher.fisherdisplacementx(rhotpt,n)
   println(io,i*dt," ",abs2(ov1)," ",abs2(ov2))
   println(io2,i*dt," ",real(qfi)," ",real(qfi2)," ",real(qfi3))
   println(io3,i*dt," ",real(neg))
 end
 end
 end
 end
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/fidelities_adiabatic.dat to see the fidelities       -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains F(t) from t=0 to ",tmax," in steps of ",tau," time units             ")
 println("   - Survival probavility                                                                           ")
 println("   - Fidelity between the state and the target state                                                ")   
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/qfi_adiabatic.dat to see the QFI                     -------------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains QFI(rho(t),A) from t=0 to ",tmax," in steps of ",tau," time units    ")
 println("   - Phase sensitivity QFI(rho(f),n)                                                                ")
 println("   - Displacement sensitivity QFI(rho(f),p)                                                         ")
 println("   - Displacement sensitivity QFI(rho(f),x)                                                         ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("--------------------------------------------------------------------------------------------------- ")
 println("-------------   Go to file output/negativity_adiabatic.dat to see the negativity volume  -----------")
 println("----------- Dynamics during the ramp                                                           -----")
 println("             The file contains delta(t) from t=0 to ",tmax," in steps of ",tau," time units         ")  
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

 
end
