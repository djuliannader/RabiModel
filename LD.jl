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
 eta = 18.0
 gamma =1.0
 psi = -1.5
 tmax = 10.0
 nn = 1000
 kk = 24
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
println("not done yet")
mswflo = wigner_eig.wigner_driven2(N,om,r,omega,gamma,eta,psi,nn,kk,L)
println("done")


ics=[0.0,0.0,2.0,0.0]
tmax=100.0


# calcula probabilidad de supervivencia
floquet=troterization.troter2(N,nn,r,om,gamma,omega,eta,psi)
cs0=dynamics.initialcoherent(ics[3],ics[4],ics[1],ics[2],hbar,N)
mensaje2 = dynamics.survivalpt(cs0,floquet,tmax,gamma)


# calcula densidad de estados 
#msj = stat.dos_rmp(N,om,r,gamma,omega,eta,psi,nn)


