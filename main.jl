push!(LOAD_PATH, pwd())
import diagonalization
import reading
import dynamics
import troterization
import stat
import wigner_eig
using LinearAlgebra



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
 nu = parse(Float64, K30)
 K31=readline(f)
 K32=readline(f)
 nn = parse(Int64, K32)
 K33=readline(f)
 K34=readline(f)
 kk = parse(Int64, K34)
 K35=readline(f)
 K36=readline(f)
 L = parse(Float64, K36)
 


##-------------
#
#
##------ converting string to a list
 lmm  = reading.stringtofloatlist(K18)
 csc0 = reading.stringtofloatlist(K24)


# printing information 
println("------------------------------------------------")
println("            Input                               ")
println("Size of the Fock space N:     ",N)
println("bosonic frequency omega:      ",om)
println("fermionic frequency R:        ",r)
println("hbar:                         ",1.0)
println("parameter delta:              ",delta)
println("carrier  parameter lambda:    ",lambda)
println("coupling parameter g  :       ",eta)
println("phase parameter psi:          ",psi)
println("frequency of the driven term: ",nu)
println("strength of the driven term : ",chi)
println("State of interest           : ",kk)
println("number of subperiods of the drive : ",nn)
println("------------------------------------------------")

# normalization of parameters
chi=chi*r

println("---------------------------------------------------------------------")
println("                             Output                                  ")
println("---------------------------------------------------------------------")
spcj=trunc(Int,(lmm[2]-lmm[1])/intl)
listaj=[lmm[1]+i*intl for i in 0:spcj]
listanu = [abs(listaj[i]) for i in 1:length(listaj)]
if flag1==1  # Spectrum
  println("------------------------------------------------")
  open("spectrum_output_g.dat","w") do file
    for j in 1:length(listaj)
      print(file,string(listaj[j]))
      evalvec=diagonalization.diagonalize(N,om,r,lambda,delta,listaj[j],psi)
      for i in evalvec[1]
       print(file,"  "," ",string(i))
      end
      println(file," ")
    end
  end
  println("See output file spectrum_output_g.dat for spectrum as a function of g")
  open("spectrum_output_lambda.dat","w") do file
  for j in 1:length(listaj)
      print(file,string(listaj[j]))
      evalvec=diagonalization.diagonalize(N,om,r,listaj[j],delta,eta,psi)
      for i in evalvec[1]
       print(file,"  "," ",string(i))
      end
      println(file," ")
  end
  end    
  println("See output file spectrum_output_lambda.dat for spectrum as a function of lambda")
  open("spectrum_quasienergies_output.dat","w") do file
  for nuj in listanu
    print(file,string(nuj))
    qs = stat.quasienergies(N,om,r,lambda,delta,nn,nuj,chi,eta,psi)
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
  for nui in listanu
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nui,eta,psi)
    vecordered = stat.orderinfvec(N,om,r,lambda,delta,nn,nui,chi,eta,psi)
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
  println("See output file spectrum_contributions_output.dat")
end



if flag1==3
   println("------ Results of the time-independent AQRM   --------------")
   message = stat.analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   evalst  = diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   wigeig  = wigner_eig.wigner_eigenstate(N,om,r,lambda,delta,eta,psi,kk,L)
   println("Ground state energy of the AQRM :",evalst[1][1])
   println("Energy of the ",kk," state of the AQRM :",evalst[1][kk])
   println("Negativity of the ",kk," state of the AQRM :",real(wigeig[1]))
   println("Purity of the ",kk," state of the AQRM     :",real(wigeig[2]))
   println("------ Results of the modulated AQRM-------------------------")
   wfloquet   = wigner_eig.wigner_driven(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L)
   rpar=stat.parameter_r(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   # mswflod  = wigner_eig.wigner_drivenqs(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L)
   println("parameter <r> of the Floquet operator : ",rpar)
   println("expectation value of the AQRM Hamiltonian in the Floquet ",kk,"-th stationary state: ",real(wfloquet[1]))
   println("Negativity of the ",kk,"-th stationary state: ",real(wfloquet[2]))
   println("Purity  of the ",kk,"-th stationary state: ",real(wfloquet[3]))
   #println("Wehlr entropy  of the ",kk,"-th stationary state: ",real(wfloquet[4]))
   #println("Maximal fidelity |<F_k|E_n>|^2= ",real(wfloquet[4])," for k= ",kk," and n= ",floor(Int,real(wfloquet[5])))
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
    mensaje4=dynamics.fotoc(cs0,tmax,1.0,N,om,r,lambda,delta,eta,psi)
    wpsit = wigner_eig.wigner_evolt(N,om,r,lambda,delta,eta,psi,L,cs0,tmax)
    println("Negativity of the state (at time=",tmax,"): ",wpsit[1])
  end
  if flag2==2 || flag2==3 # Floquet
    println("----   Dynamics of the modulated AQRM initiates     ----")
    tau = 2*pi/nu
    pf = floor(tmax/tau)  
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
    mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,nu)
    mensaje3 = dynamics.fotoct(cs0,floquet,tmax,nu,N)
    wpsit_floquet = wigner_eig.wigner_evolt_driven(N,om,r,lambda,delta,eta,psi,nu,chi,nn,L,cs0,pf)
    println("Negativity of the state (at time=",pf*tau,"): ",wpsit_floquet[1])
  end
end



end