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
  pop = (im/(2)^(1/2))*( adop -  aop)
  #pop = (1/(2)^(1/2))*( adop +  aop)
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


function fishern(rho,Nmax)
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
       sd = lstatead*nop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst = 2*((lam[l]-lam[k])^2/(lam[k]+lam[l]))*abs2(sd) 
         qfi = qfi + qfiinst
	 #println("here ",qfiinst)
       end
    end
  end
  bc=FockBasis(Nmax)
  nn = number(bc) 
  nexp = real(expect(rho,nn))
  return real(qfi)/(4*nexp)
end

function fishern2(rho,Nmax)
  rhodiag = eigenstates(dense(rho))
  nop = numberop(Nmax)
  lam = rhodiag[1]
  qfi2t = 0 + 0*im
  qfi1t = 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:k-1
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*nop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])^2/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     for l in k+1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*nop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])^2/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
     kstatead = transpose(conj(kstate))
     n2ev = kstatead*(nop^2)*kstate
     nev2 = abs2(kstatead*(nop)*kstate)
     qfiinst1t = 4*lam[k]*(n2ev -nev2)
     qfi1t = qfi1t + qfiinst1t
  end
  bc=FockBasis(Nmax)
  nn = number(bc) 
  nexp = real(expect(rho,nn))
  return real(qfi1t - qfi2t)/(4*(nexp))
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