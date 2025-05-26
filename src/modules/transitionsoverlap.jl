module transitionoverlap
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
import DQPT




n=110
om=1.0
ri= 0.0
rf= 400.0
rint = 1.0
gi =  2.0
gf =  0.5
gint = 0.1
lambda0 = 3/20
delta0 = 0.0
psi0 =0.0
alpha=1.0
ph=0.0


name = "densityoverlap.dat"

igmax = convert(Int,(gi-gf)/gint)
irmax = convert(Int,(rf-ri)/rint)


open(name,"w") do io
rinst = ri
for ir in 1:irmax
  ginst = gi
  for ig in 1:igmax
    istate = DQPT.initialstatequench(n,om,rinst,lambda0,delta0,gi,psi0)
    phi0 = alpha^(1/2)*istate[1] + (1-alpha)^(1/2)*exp(im*ph)*istate[2]
    hamf = diagonalization.hamiltonian(n,om,rinst,lambda0,delta0,ginst,psi0)
    ovl = DQPT.overlapdqptlist(phi0,hamf,n)
    k0 = ovl[1]
    kmax = maximum(ovl)
    println(io,ginst-gi," ",rinst," ",k0/kmax)
    ginst = ginst - gint
  end
  rinst = rinst + rint
end

end



end
