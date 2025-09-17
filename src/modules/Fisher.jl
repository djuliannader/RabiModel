push!(LOAD_PATH, pwd())
module Fisher
using LinearAlgebra
using QuantumOptics
using HCubature  # loaded by CalculusWithJulia
using PyPlot


function fisherdisplacementp(rho,Nmax)
  rhodiag = eigenstates(dense(rho))
  aop  = anhilation(Nmax)
  adop = transpose(aop)
  sqop = (1/(2^(1/2)))*(aop+adop)
  lam = rhodiag[1]
  qfi2t = 0 + 0*im
  qfi1t = 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:k-1
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     for l in k+1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
     kstatead = transpose(conj(kstate))
     sq2ev = kstatead*(sqop^2)*kstate
     sqev2 = abs2(kstatead*(sqop)*kstate)
     qfiinst1t = 4*lam[k]*(sq2ev -sqev2)
     qfi1t = qfi1t + qfiinst1t
  end
  bc=FockBasis(Nmax)
  nn = number(bc) 
  nexp = real(expect(rho,nn))
  return real(qfi1t - qfi2t)/2
  #return real(qfi1t - qfi2t)/(2*(1+4*nexp))
end

function fisherdisplacementx(rho,Nmax)
  rhodiag = eigenstates(dense(rho))
  aop  = anhilation(Nmax)
  adop = transpose(aop)
  sqop = (1/(im*2^(1/2)))*(aop-adop)
  lam = rhodiag[1]
  qfi2t = 0 + 0*im
  qfi1t = 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:k-1
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     for l in k+1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
     kstatead = transpose(conj(kstate))
     sq2ev = kstatead*(sqop^2)*kstate
     sqev2 = abs2(kstatead*(sqop)*kstate)
     qfiinst1t = 4*lam[k]*(sq2ev -sqev2)
     qfi1t = qfi1t + qfiinst1t
  end
  bc=FockBasis(Nmax)
  nn = number(bc) 
  nexp = real(expect(rho,nn))
  return real(qfi1t - qfi2t)/2
end

function fishersqueezing(rho,Nmax)
  rhodiag = eigenstates(dense(rho))
  sqop = squeezingop(Nmax)
  lam = rhodiag[1]
  qfi2t = 0 + 0*im
  qfi1t = 0 + 0*im
  epsilon=0.0000001
  for k in 1:length(lam)
     for l in 1:k-1
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     for l in k+1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*sqop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
     kstatead = transpose(conj(kstate))
     sq2ev = kstatead*(sqop^2)*kstate
     sqev2 = abs2(kstatead*(sqop)*kstate)
     qfiinst1t = 4*lam[k]*(sq2ev -sqev2)
     qfi1t = qfi1t + qfiinst1t
  end
  return real(qfi1t - qfi2t)
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
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     for l in k+1:length(lam)
       kstate = [rhodiag[2][k][i] for i in 1:length(lam)]
       lstate = [rhodiag[2][l][i] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*nop*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
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
  #println("------->here:",nexp)
  return real(qfi1t - qfi2t)/(4*(nexp))
end


function anhilation(Nmax)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end

function xoperator(Nmax)
  aop  = anhilation(Nmax)
  adop = transpose(aop)
  xop = (1/2^(1/2))*(aop+adop)
  return xop
end

function numberop(Nmax)
  diag  = [(i-1) for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [0.0 for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end

function squeezingop(Nmax)
  aop=anhilation(Nmax)
  adop = transpose(aop)
  sqop = adop^2-aop^2
  return sqop
end


end
