push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/diagonalization.jl")
include("modules/reading.jl")
include("modules/stat.jl")
include("modules/wigner_eig.jl")
include("modules/dynamics.jl")
include("modules/troterization.jl")
using .diagonalization
using .troterization
using .stat
using .wigner_eig
using .dynamics



println("\r Rabi Model ")
println("\r Initiating ")




# Reading data from input file
#------------------------------------
open("input.dat") do f
 K1=readline(f)
 K2=readline(f)
 N = parse(Int64, K2)
 K3=readline(f)
 K4=readline(f)
 om = parse(Float64, K4)
 K5=readline(f)
 K6=readline(f)
 r = parse(Float64, K6)
 K7=readline(f)
 K8=readline(f)
 delta = parse(Float64, K8)
 K9=readline(f)
 K10=readline(f)
 lambda = parse(Float64, K10)
 K11=readline(f)
 K12=readline(f)
 eta = parse(Float64, K12)
 K13=readline(f)
 K14=readline(f)
 psi = parse(Float64, K14)
 K15=readline(f)
 K16=readline(f)
 flag1 = parse(Int64, K16)
 K17=readline(f)
 K18=readline(f)
 K19=readline(f)
 K20=readline(f)
 intl = parse(Float64, K20)
 K21=readline(f)
 K22=readline(f)
 flag2  = parse(Int64, K22)
 K23=readline(f)
 K24=readline(f)
 K25=readline(f)
 K26=readline(f)
 tmax = parse(Float64, K26)
 K27=readline(f)
 K28=readline(f)
 chi = parse(Float64, K28)
 K29=readline(f)
 K30=readline(f)
 tau = parse(Float64, K30)
 K31=readline(f)
 K32=readline(f)
 nn = parse(Int64, K32)
 K33=readline(f)
 K34=readline(f)
 kk = parse(Int64, K34)
 K35=readline(f)
 K36=readline(f)
 L = parse(Float64, K36)
 K37=readline(f)
 K38=readline(f)
 flagt = parse(Int64, K38)
 K39=readline(f)
 K40=readline(f)
 lk = parse(Int64, K40)

##-------------
#
#
##------ converting string to a list
 lmm  = reading.stringtofloatlist(K18)
 csc0 = reading.stringtofloatlist(K24)



# printing information 
println("------------------------------------------------")
println("            Input                               ")
println("Size of the Fock space   (N)         : ",N)
println("Oscillator frequency   (Omega)       : ",om)
println("Qubit frequency          (R)         : ",r)
println("QRM (0), JC (1) or AJC (-1)          : ",delta)
println("Carrier  parameter     (lambda)      : ",lambda)
println("Coupling parameter       (g)         : ",eta)
println("Phase parameter         (psi)        : ",psi)
println("Modulation period       (Tau)        : ",tau)
println("Modulation amplitude     (Xi)        : ",chi)
println("State of interest        (k)         : ",kk)
println("Number of subperiods of the modulation    : ",nn)
println("------------------------------------------------")


# normalization of parameters
chi=chi*r
nu=2*pi/tau

println("---------------------------------------------------------------------")
println("                             Output                                  ")
println("---------------------------------------------------------------------")
spcj=trunc(Int,(lmm[2]-lmm[1])/intl)
listaj=[lmm[1]+i*intl for i in 0:spcj]
listanu = [abs(listaj[i]) for i in 1:length(listaj)]
if flag1==1  # Spectrum
  println("------------------------------------------------")
  open("output/spectrum_output_g.dat","w") do file
    for j in 1:length(listaj)
      print(file,string(listaj[j]))
      evalvec=diagonalization.diagonalize(N,om,r,lambda,delta,listaj[j],psi)
      for i in evalvec[1]
       print(file,"  "," ",string(i))
      end
      println(file," ")
    end
  end
  println("See output file output/spectrum_output_g.dat for spectrum as a function of g")
  open("output/spectrum_output_lambda.dat","w") do file
  for j in 1:length(listaj)
      print(file,string(listaj[j]))
      evalvec=diagonalization.diagonalize(N,om,r,listaj[j],delta,eta,psi)
      for i in evalvec[1]
       print(file,"  "," ",string(i))
      end
      println(file," ")
  end
  end
  println("See output file output/spectrum_output_lambda.dat for spectrum as a function of lambda")
  if lmm[1]>0.0
  open("output/spectrum_quasienergies_output.dat","w") do file
  for nuj in listanu
    print(file,string(nuj))
    qs = stat.quasienergies(N,om,r,lambda,delta,nn,nuj,chi,eta,psi,flagt)
    for i in qs
     print(file,"  "," ",string(i))
    end
    println(file," ")
  end
  end
  println("See output file spectrum_quasienergies_output.dat")
   h0=diagonalization.hamiltonian(N,om,r,lambda,delta,eta,psi)
   eigvecsh0 = eigvecs(h0)
   open("spectrum_contributions_output.dat","w") do io
   ovlist=[0.0,0.0,0.0,0.0,0.0]
   #println("------------>here")
   for nui in listanu
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nui,eta,psi,flagt)
    vecordered = stat.orderinfvec(N,om,r,lambda,delta,nn,nui,chi,eta,psi,flagt)
    evsf=eigvecs(floquet)
    listfstatek = [evsf[j,vecordered[kk]] for j in 1:2*N]
    psifk = wigner_eig.buildingstate(listfstatek,N)
    for kh0 in 1:5
      listh0statek = [eigvecsh0[i,kh0] for i in 1:2*N]
      lhk = transpose(listh0statek)
      ov = abs2(lhk*listfstatek)
      ovlist[kh0]=ov
    end
    println(io,nui,"  ",ovlist[1]," ",ovlist[2]," ",ovlist[3]," ", ovlist[4]," ",ovlist[5])
  end
  end
  println("See output file output/spectrum_contributions_output.dat")
  end
end



if flag1==3
   println("------ Results of the time-independent AQRM   --------------")
   message = stat.analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   evalst  = diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   wigeig  = wigner_eig.wigner_eigenstate(N,om,r,lambda,delta,eta,psi,kk,L)
   println("Ground state energy of the AQRM                  : ",evalst[1][1])
   println("Expectation value <E_k|H_0|E_k>                  : ",evalst[1][kk])
   println("Bosonic Negativity of the ",kk," state of the AQRM   delta(rho_b)    : ",real(wigeig[1]))
   println("Purity of the ",kk," state of the AQRM           : ",real(wigeig[2]))
    println("------ Results of AQRM with a small modulation in the qubit frequency ------")
    println("---------------------------------------------------------------------------")
   wfloquet   = wigner_eig.wigner_driven(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L,flagt,lk)
   rpar=stat.parameter_r(N,om,r,lambda,delta,nn,nu,chi,eta,psi,flagt)
   # mswflod  = wigner_eig.wigner_drivenqs(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L)
   println("Parameter <r> of the Floquet operator            : ",rpar)
   println("Expectation value      <F_k|H_0|F_k>             : ",real(wfloquet[1]))
   println("Bosonic Negativity of the ",kk,"-th stationary state  delta(rho_b)  : ",real(wfloquet[2]))
   println("Purity  of the ",kk,"-th stationary state        : ",real(wfloquet[3]))
   #println("Wehlr entropy  of the ",kk,"-th stationary state: ",real(wfloquet[4]))
   #println("Maximal fidelity |<F_k|E_n>|^2= ",real(wfloquet[4])," for k= ",kk," and n= ",floor(Int,real(wfloquet[5])))
   for jk in 1:lk
     println("Overlap |<E_l|F_k>|^2 for k= ",kk," and l=",jk,"  : ",real(wfloquet[4][jk]))
   end
   for jk in 1:lk
     println("Overlap |<n|F_k>|^2 for k= ",kk," and n=",jk,"  : ",real(wfloquet[5][jk]))
   end
   println("See file levels_output.dat")
   println("-------------------------------------------------------------")
   # --------------------------------------
  end


if flag1==2  # Dynamics
  cs0=dynamics.initialcoherent(csc0[3],csc0[4],csc0[1],csc0[2],1.0,N)
  #cs0=dynamics.initialsqueezed(csc0[3],csc0[4],csc0[1],csc0[2],1.0,N,0.5,pi/2)
  if flag2==1 || flag2==3 # static Hamiltonian
    println("---- Dynamics of the time independent AQRM initiates ----")
    mensaje1=dynamics.survivalp(cs0,tmax,1.0,N,om,r,lambda,delta,eta,psi)
    av=dynamics.fotoc(cs0,tmax,1.0,N,om,r,lambda,delta,eta,psi,L)
    #wpsit = wigner_eig.wigner_evolt(N,om,r,lambda,delta,eta,psi,L,cs0,tmax)
    #println("Negativity of the state (at time=",tmax,"): ",wpsit[1])
    println("------------ Averages ---------------")
    println("Average Negativity volume over the period T=",tmax," : ",av[1])
    println("Average QFI[rho(t),n]/(4n) over the period T=",tmax," : ",av[2])
    println("Average QFI[rho(t),p]/2 over the period T=",tmax," : ",av[3])
    println("Average QFI[rho(t),x]/2 over the period T=",tmax," : ",av[4])
    println("Average QFI[rho,xp+px]/(4n+2) over the period T=",tmax," : ",av[5])
    if flag2==1  
      println("---------------------------------------")
      println(" Wigner function at time t=",tmax)
    end
  end
  if flag2==2 || flag2==3 # Floquet
    println("----   Dynamics of the modulated AQRM initiates     ----")
    pf = floor(tmax/tau)  
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi,flagt)
    mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,nu)
    av = dynamics.fotoct(cs0,floquet,tmax,nu,r,N,L)
    #wpsit_floquet = wigner_eig.wigner_evolt_driven(N,om,r,lambda,delta,eta,psi,nu,chi,nn,L,cs0,pf,flagt)
    #println("Negativity of the state (at time=",pf*tau,"): ",wpsit_floquet[1])
    println("------------ Averages ---------------")
    println("Average Negativity volume over the period T=",tmax," : ",av[1])
    println("Average QFI[rho(t),n]/(4n) over the period T=",tmax," : ",av[2])
    println("Average QFI[rho(t),p]/2 over the period T=",tmax," : ",av[3])
    println("Average QFI[rho(t),x]/2 over the period T=",tmax," : ",av[4])
    println("Average QFI[rho,xp+px]/(4n+2) over the period T=",tmax," : ",av[5])
    println("---------------------------------------")
    println(" Wigner function at time t=",pf*tau)
  end 
  psit = av[6]
  rhot = psit*transpose(conj(psit))  
  wig = wigner_eig.wigner_rhot(rhot,L,r,N)
  println("---------------------------------------")
end



end
