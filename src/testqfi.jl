module testqfi
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
import Fisher


Nmax=100
b = FockBasis(Nmax)
al1 = (1/2^(1/2))*(4 + 0*im)
cs1 = coherentstate(b,al1)
al2 = (1/2^(1/2))*(-4 + 0*im)
cs2 = coherentstate(b,al2)
rho1 = dm(cs1)
rho2 = dm(cs2)
#rho = (1/2)*rho1 + (1/2)*rho2
cat = (1/2^(1/2))*(cs1+cs2)
rho=dm(cat)
qfi  = Fisher.fisherdisplacementp(rho,Nmax)

println("qfi :",qfi)

end