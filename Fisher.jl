push!(LOAD_PATH, pwd())
module Fisher
using LinearAlgebra
using QuantumOptics
using HCubature  # loaded by CalculusWithJulia
using PyPlot
export fisherp



function fisherp(rho,Nmax)
  rhodiag = eigenstates(dense(rho))
  aop = creation(Nmax)
  adop = transpose(aop)
  pop = (im/(2)^(1/2))*(adop - aop)
  lam = rhodiag[1]
  qfi= 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*pop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst = 2*((lam[l]-lam[k])^2/(lam[k]+lam[l]))*abs2(sd) 
         qfi = qfi + qfiinst
	 #println("here ",qfiinst)
       end
    end
  end
  return real(qfi)
end


function fishern(rho,Nmax,phi)
  rhodiag = eigenstates(dense(rho))
  nop = numberop(Nmax)
  lam = rhodiag[1]
  qfi= 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*exp(im*nop*phi)*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst = 2*((lam[l]-lam[k])^2/(lam[k]+lam[l]))*abs2(sd) 
         qfi = qfi + qfiinst
	 #println("here ",qfiinst)
       end
    end
  end
  return real(qfi)
end

function fishervsphi(rho,Nmax)
   philist = [i*0.1 for i in 0:62]
   #println("flag:",length(philist))
   open("qfivsphi.dat","w") do io
   for i in philist
     println(i/(2*pi))
     qfi = fishern(rho,Nmax,i)
     println(io,i," ",qfi)
   end
   println(" --   Go to file qfivsphi.dat to see the QFI[rho,exp(i*n*phi)] as a function of phi  --")
   end
end


function creation(Nmax)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end

function numberop(Nmax)
  diag  = [(i-1) for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [0.0 for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end


end