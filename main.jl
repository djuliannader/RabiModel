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
 hbar = parse(Float64, K8)
 K9=readline(f)
 K10=readline(f)
 delta = parse(Float64, K10)
 K11=readline(f)
 K12=readline(f)
 lambda = parse(Float64, K12)
 K13=readline(f)
 K14=readline(f)
 eta = parse(Float64, K14)
 K15=readline(f)
 K16=readline(f)
 psi = parse(Float64, K16)
 K17=readline(f)
 K18=readline(f)
 flag1 = parse(Int64, K18)
 K19=readline(f)
 K20=readline(f)
 K21=readline(f)
 K22=readline(f)
 intl = parse(Float64, K22)
 K23=readline(f)
 K24=readline(f)
 flag2  = parse(Int64, K24)
 K25=readline(f)
 K26=readline(f)
 K27=readline(f)
 K28=readline(f)
 tmax = parse(Float64, K28)
 K29=readline(f)
 K30=readline(f)
 chi = parse(Float64, K30)
 K31=readline(f)
 K32=readline(f)
 nu = parse(Float64, K32)
 K33=readline(f)
 K34=readline(f)
 nn = parse(Int64, K34)
 K35=readline(f)
 K36=readline(f)
 kk = parse(Int64, K36)
 K37=readline(f)
 K38=readline(f)
 L = parse(Int64, K38)
 

##-------------
#
#
##------ converting string to a list
 lmm  = reading.stringtofloatlist(K20)
 csc0 = reading.stringtofloatlist(K26)

# printing information 
println("------------------------------------------------")
println("            Input                               ")
println("Size of the Fock space N:     ",N)
println("bosonic frequency omega:      ",om)
println("fermionic frequency R:        ",r)
println("hbar:                         ",hbar)
println("parameter delta:              ",delta)
println("carrier  parameter lambda:    ",lambda)
println("coupling parameter eta:       ",eta)
println("phase parameter psi:          ",psi)
println("frequency of the driven term: ",nu)
println("strength of the driven term : ",chi)
println("number of subperiods of the driving : ",nn)
println("------------------------------------------------")


println("------------------------------------------------")
println("            Output                              ")
if flag1==1  # Spectrum
  println("See output file spectrum_output.dat")
  println("------------------------------------------------")
  open("spectrum_output.dat","w") do file
  spcj=trunc(Int,(lmm[2]-lmm[1])/intl)
  listaj=[lmm[1]+i*intl for i in 0:spcj]
  for j in listaj
    print(file,string(j))
    evalvec=diagonalization.diagonalize(N,om,r,j,delta,eta,psi)
    for i in evalvec[1]
      print(file,"  "," ",string(i))
    end
  println(file," ")
  end
  end
  nup=1.0  # initial frequency
  nuint=0.02
  open("spectrum_quasienergies_output.dat","w") do file
  for i in 1:50
  print(file,string(nup))
  qs = stat.quasienergies(N,om,r,lambda,delta,nn,nup,chi,eta,psi)
  for i in qs
      print(file,"  "," ",string(i))
    end
    println(file," ")
  nup=nup+nuint
  end
  end
  end


if flag1==3
   message=stat.analysisH(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   rpar=stat.parameter_r(N,om,r,lambda,delta,nn,nu,chi,eta,psi)
   evalst=diagonalization.diagonalize(N,om,r,lambda,delta,eta,psi)
   mswigeig = wigner_eig.wigner_eigenstate(N,om,r,lambda,delta,eta,psi,kk,L)
   mswflo   = wigner_eig.wigner_driven(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L)
   # mswflod  = wigner_eig.wigner_drivenqs(N,om,r,lambda,delta,eta,psi,nu,chi,nn,kk,L)
   println("Ground state energy of the time independent Hamiltonian :",evalst[1][1])
   println("parameter <r> of the Floquet operator :",rpar)
   println("See file levels_output.dat")
   # --------------------------------------
  end


if flag1==2  # Dynamics
  cs0=dynamics.initialcoherent(csc0[3],csc0[4],csc0[1],csc0[2],hbar,N)
  if flag2==1 || flag2==3 # static Hamiltonian
    mensaje1=dynamics.survivalp(cs0,tmax,hbar,N,om,r,lambda,delta,eta,psi)
  end
  if flag2==2 || flag2==3 # Floquet
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu,eta,psi)
    mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,nu)
  end
end



end