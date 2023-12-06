push!(LOAD_PATH, pwd())
import diagonalization
import reading
import dynamics

println("\r Initiating spectrum ")

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
 flag1 = parse(Int64, K12)
 K13=readline(f)
 K14=readline(f)
 K15=readline(f)
 K16=readline(f)
 intl = parse(Float64, K16)
 K17=readline(f)
 K18=readline(f)
 lambda = parse(Float64, K18)
 K19=readline(f)
 K20=readline(f)
 K21=readline(f)
 K22=readline(f)
 tmax = parse(Float64, K22)

#-------------


#------ converting string to a list
 lmm  = reading.stringtofloatlist(K14)
 csc0 = reading.stringtofloatlist(K20)

# printing information 
println("------------------------------------------------")
println("Size of the Fock space N=",N)
println("parameter omega=",om)
println("parameter=",r)



if flag1==1  # Spectrum
  println("See output file spectrum_output.dat")
  println("------------------------------------------------")
  open("spectrum_output.dat","w") do file
  spcj=trunc(Int,(lmm[2]-lmm[1])/intl)
  listaj=[lmm[1]+i*intl for i in 0:spcj]
  for j in listaj
    print(file,string(j))
    evalvec=diagonalization.diagonalize(50,1,0.8,j,1)
    for i in evalvec[1]
      print(file,"  "," ",string(i))
    end
  println(file," ")
  end
  end
  end


if flag1==2  # Dynamics
  cs0=dynamics.initialcoherent(csc0[3],csc0[4],csc0[1],csc0[2],hbar,N)
  mensaje=dynamics.survivalp(cs0,tmax,hbar,N,om,r,lambda,delta)
end


end