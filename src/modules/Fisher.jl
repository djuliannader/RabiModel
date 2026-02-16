module Fisher
push!(LOAD_PATH, pwd())
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


function fisherdisplacement(rho,Nmax,phi)
  rhodiag = eigenstates(dense(rho))
  aop  = anhilation(Nmax)
  adop = transpose(aop)
  sqop = (1/(2^(1/2)))*(adop*exp(im*phi)+aop*exp(-im*phi))
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
  bc=FockBasis(Nmax)
  nn = number(bc) 
  nexp = real(expect(rho,nn))
  #println("------->here:",nexp)
  return real(qfi1t - qfi2t)/(4*(nexp+1/2))
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
  sqop = (im/2)*(adop^2-aop^2)
  return sqop
end


function cfisherdisplacement2(rhoqo,Nmax,phi)
  theta= 0.3
  delta=0.1  
  rhodiag = eigenstates(dense(rhoqo))
  fb = FockBasis(Nmax)    
# -------------------------------
  xmin = 5.0          # should be Float64
  Nm   = 200
  bx = PositionBasis(-xmin, xmin, Nm)
  x  = samplepoints(bx)   # x_i points
  dx = spacing(bx)        # Δx  
  Txf  = transform(bx, fb)
    # --------------
  aop = destroy(fb)
  adop = dagger(aop)  
  gop = (1/2^(1/2))*(adop*exp(im*phi)+aop*exp(-im*phi))
  nop = adop*aop  
  lam = rhodiag[1]
  global cfi = 0.0 + 0.0*im
  rhot2 = exp(im*gop*(theta+delta))*rhoqo*exp(-im*gop*(theta+delta))
  rhot1 = exp(im*gop*(theta-delta))*rhoqo*exp(-im*gop*(theta-delta))
  rhoqor = exp(im*gop*(theta))*rhoqo*exp(-im*gop*(theta))    
  rhot2p = exp(im*nop*(phi+pi/2))*rhot2*exp(-im*nop*(phi+pi/2))
  rhot1p = exp(im*nop*(phi+pi/2))*rhot1*exp(-im*nop*(phi+pi/2))  
  rhoqorp = exp(im*nop*(phi+pi/2))*rhoqor*exp(-im*nop*(phi+pi/2))    
  rhot2x = dense(Txf * rhot2p * dagger(Txf))
  rhot1x = dense(Txf * rhot1p * dagger(Txf))  
  rhox = dense(Txf * rhoqorp * dagger(Txf))  
  qt2 = real.(diag(rhot2x.data))   # probability mass
  pt2 = qt2 ./ dx
  qt1 = real.(diag(rhot1x.data))   # probability mass
  pt1 = qt1 ./ dx  
  q = real.(diag(rhox.data))   # probability mass
  p = q ./ dx  
  for i in 1:length(x)
     global cfi = cfi + ((pt2[i]-pt1[i])/(2*delta))^2*dx/p[i]
  end      
  return real(cfi)/2
end


function cfisherdisplacement(rhoqo,Nmax,theta)
  phi= theta + pi/2
  fb = FockBasis(Nmax)    
# -------------------------------
  aop = destroy(fb)
  adop = dagger(aop)  
  Gop = (1/2^(1/2))*(aop*exp(im*theta) + adop*exp(-im*theta))
  Mop = (1/2^(1/2))*(aop*exp(im*phi) + adop*exp(-im*phi))
  expn = expect(rhoqo,Gop*Mop - Mop*Gop)
  num = abs2(-im*expn)  
  M2exp = expect(rhoqo,Mop*Mop)  
  Mexp = expect(rhoqo,Mop)
  VarM = M2exp - Mexp^2
  cfi = num/VarM  
  return real(cfi)/2  
end

function cfisherrotationb(rhoqo,Nmax,theta,beta)
    fb = FockBasis(Nmax)
    aop = destroy(fb)
    adop = dagger(aop)  
    nop = adop*aop
    rhodiag = eigenstates(dense(rhoqo))
    lam = rhodiag[1]
    rhotheta = exp(im*theta*nop)*rhoqo*exp(-im*theta*nop)
    rhobeta = exp(beta*adop - conj(beta)*aop)*rhotheta*exp(conj(beta)*aop - beta*adop)
    rhobdat = rhobeta.data
    drho = im*(nop*rhotheta - rhotheta*nop)
    drhob = exp(beta*adop - conj(beta)*aop)*drho*exp(conj(beta)*aop - beta*adop)
    drhobdat = drhob.data
    sum = 0.0 + 0*im
    for n in 1:length(lam)
        num = abs2(drhobdat[n,n])
        den = rhobdat[n,n]
       sum  = sum + num/den
    end
    nexp = real(expect(rhoqo,nop))
    return real(sum)/(4*nexp)
end

#Nmax=50
#fb = FockBasis(Nmax)
#println("here")
#psi0 = fockstate(fb,0)
#al = 2
#psi1 = coherentstate(fb,al)
#psi2 = coherentstate(fb,-al)
#psi0 = (1/(2*(1+exp(-2*abs2(al))))^(1/2))*(psi1+psi2)
#rho0 = dm(psi0)
#beta=[i for i in 0:0.5:20]
#open("test.dat","w") do io
#for b in beta
#    cc = cfisherdisplacement2(rho0,Nmax,pi/8)
#    println("result: ",cc)
#end
#end


# Displacement photon counting
function cfisherrotation(rhoqo,Nmax)
    beta =  [i + j*im for i in 0.1:0.1:5 for j in 0.1:0.1:5]
    cfi=[]
    for b in beta
        cc = cfisherrotationb(rhoqo,Nmax,0.5,b)
        append!(cfi,cc)
    end
    val, idx = findmax(cfi)
    return [beta[idx],val]
end

end
