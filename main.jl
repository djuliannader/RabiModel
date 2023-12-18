push!(LOAD_PATH, pwd())
import diagonalization
import reading
import dynamics
import troterization
import statistics

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
 flag1 = parse(Int64, K14)
 K15=readline(f)
 K16=readline(f)
 K17=readline(f)
 K18=readline(f)
 intl = parse(Float64, K18)
 K19=readline(f)
 K20=readline(f)
 flag2  = parse(Int64, K20)
 K21=readline(f)
 K22=readline(f)
 K23=readline(f)
 K24=readline(f)
 tmax = parse(Float64, K24)
 K25=readline(f)
 K26=readline(f)
 chi = parse(Float64, K26)
 K27=readline(f)
 K28=readline(f)
 nu = parse(Float64, K28)
 K29=readline(f)
 K30=readline(f)
 nn = parse(Int64, K30)
 

#-------------


#------ converting string to a list
 lmm  = reading.stringtofloatlist(K16)
 csc0 = reading.stringtofloatlist(K22)

# printing information 
println("------------------------------------------------")
println("Size of the Fock space N:  ",N)
println("bosonic frequency omega:   ",om)
println("fermionic frequency R:     ",r)
println("hbar:                      ",hbar)
println("coupling parameter lambda: ",lambda)
println("parameter delta:           ",delta)



if flag1==1  # Spectrum
  println("See output file spectrum_output.dat")
  println("------------------------------------------------")
  open("spectrum_output.dat","w") do file
  spcj=trunc(Int,(lmm[2]-lmm[1])/intl)
  listaj=[lmm[1]+i*intl for i in 0:spcj]
  for j in listaj
    print(file,string(j))
    evalvec=diagonalization.diagonalize(N,om,r,2*j,delta)
    for i in evalvec[1]
      print(file,"  "," ",string(i))
    end
  println(file," ")
  end
  end
  end


if flag1==3
   message=statistics.analysisH(N,om,r,lambda,delta,nn,nu,chi)
   println("See file levels_output.dat")
end


if flag1==2  # Dynamics
  cs0=dynamics.initialcoherent(csc0[3],csc0[4],csc0[1],csc0[2],hbar,N)
  if flag2==1 || flag2==3 # static Hamiltonian
    mensaje1=dynamics.survivalp(cs0,tmax,hbar,N,om,r,lambda,delta)
  end
  if flag2==2 || flag2==3 # Floquet
    floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu)
    mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,nu)
  end
end


end