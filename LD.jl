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
 N = 50
 om = 1.0
 r = 50.0
 hbar = 1.0
 omega = 1.0 
 eta = 15.0
 gamma = 0.5
 psi = -pi/2
 tmax = 10.0
 nn = 1000
 kk = 9
 L = 10
 



# printing information 
println("------------------------------------------------")
println("            Input                               ")
println("Size of the Fock space N:     ",N)
println("bosonic frequency omega:      ",om)
println("fermionic frequency R:        ",r)
println("hbar:                         ",hbar)
println("carrier  parameter Omega:    ",omega)
println("parameter gamma         :    ",gamma)
println("coupling parameter eta:       ",eta)
println("phase parameter psi:          ",psi)
println("number of subperiods of the driving : ",nn)
println("------------------------------------------------")


# calcula funcion de wigner de estados estacionarios
println("Calculating Wigner function of the ",kk," stationary state")
mswflo = wigner_eig.wigner_driven2(N,om,r,omega,gamma,eta,psi,nn,kk,L)
#mswflo = wigner_eig.wigner_driven3(N,om,r,omega,gamma,eta,psi,nn,kk,L,8)
println("done")

# calcula paridad y entropia de entrelazamiento
println("Calculating purity of the stationary state")
plist = stat.purity_rmp(N,om,r,omega,gamma,eta,psi,nn)
#lambda=0.01
#delta=0
#eta=12.0
#psi=0.0
#plist = stat.purity(N,om,100,lambda,delta,eta,psi)
#println(plist[1])
println(plist[1])
println(plist[2])


ics=[0.0,0.0,4.0,0.0]
tmax=100.0


# calcula probabilidad de supervivencia
floquet=troterization.troter2(N,nn,r,om,gamma,omega,eta,psi)
floquet2=troterization.troter3(N,nn,r,om,gamma,omega,eta,psi,4)
cs0=dynamics.initialcoherent(ics[3],ics[4],ics[1],ics[2],hbar,N)
mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,gamma)
mensaje3 = dynamics.survivalpt2(cs0,floquet2,tmax,gamma)


# calcula densidad de estados 
msj = stat.dos_rmp(N,om,r,gamma,omega,eta,psi,nn)


