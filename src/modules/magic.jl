push!(LOAD_PATH, pwd())
module magic
using LinearAlgebra
using QuantumOptics
using Convex, SCS
using PyPlot
export robustness, dicrete_wigner

function robustness(rho)

#   Spin basis
    ba = SpinBasis(1//2)

#   Identity
    idop = dm(spindown(ba)) + dm(spinup(ba))

    
#   bi
    b1 = tr(idop*rho)
    b2 = tr(sigmax(ba)*rho)
    b3 = tr(sigmay(ba)*rho)
    b4 = tr(sigmaz(ba)*rho)

    b = [real(b1),real(b2),real(b3),real(b4)]


#   qubit stabilizer state projectors    
    s1 = dm((1/2^(1/2))*(spindown(ba) + spinup(ba)))
    s2 = dm((1/2^(1/2))*(spindown(ba) - spinup(ba)))
    s3 = dm((1/2^(1/2))*(spindown(ba) +im*spinup(ba)))
    s4 = dm((1/2^(1/2))*(spindown(ba) -im*spinup(ba)))
    s5 = dm(spindown(ba))
    s6 = dm(spinup(ba))

    
#   Aji    

    a11 = tr(idop*s1)
    a21 = tr(idop*s2)
    a31 = tr(idop*s3)
    a41 = tr(idop*s4)
    a51 = tr(idop*s5)
    a61 = tr(idop*s6)

    a12 = tr(sigmax(ba)*s1)
    a22 = tr(sigmax(ba)*s2)
    a32 = tr(sigmax(ba)*s3)
    a42 = tr(sigmax(ba)*s4)
    a52 = tr(sigmax(ba)*s5)
    a62 = tr(sigmax(ba)*s6)

    a13 = tr(sigmay(ba)*s1)
    a23 = tr(sigmay(ba)*s2)
    a33 = tr(sigmay(ba)*s3)
    a43 = tr(sigmay(ba)*s4)
    a53 = tr(sigmay(ba)*s5)
    a63 = tr(sigmay(ba)*s6)

    a14 = tr(sigmaz(ba)*s1)
    a24 = tr(sigmaz(ba)*s2)
    a34 = tr(sigmaz(ba)*s3)
    a44 = tr(sigmaz(ba)*s4)
    a54 = tr(sigmaz(ba)*s5)
    a64 = tr(sigmaz(ba)*s6)
 
A=[real(a11) real(a21) real(a31) real(a41) real(a51) real(a61);
   real(a12) real(a22) real(a32) real(a42) real(a52) real(a62);
   real(a13) real(a23) real(a33) real(a43) real(a53) real(a63);
    real(a24) real(a24) real(a34) real(a44) real(a54) real(a64)]

#  Solving the nxm linear system 
   x = Variable(size(A,2))
   problem = minimize(norm(x,1), A*x == b)
   solve!(problem, SCS.Optimizer; silent = true)
   x_opt = evaluate(x);

   return sum(abs.(x_opt))
    
end

function discrete_wigner(rho)

#   Spin basis
    ba = SpinBasis(1//2)

#   Identity
    idop = dm(spindown(ba)) + dm(spinup(ba))


    idop = dm(spindown(ba)) + dm(spinup(ba)) 
    A00 = (1/2)*(idop + sigmax(ba) + sigmay(ba) + sigmaz(ba)) 
    A01 = (1/2)*(idop + sigmax(ba) - sigmay(ba) - sigmaz(ba))
    A10 = (1/2)*(idop - sigmax(ba) + sigmay(ba) - sigmaz(ba))
    A11 = (1/2)*(idop - sigmax(ba) - sigmay(ba) + sigmaz(ba)) 

    W = [0.0   0.0;
        0.0   0.0]

    W[1,1] = real((1/2)*tr(rho*A00))
    W[1,2] = real((1/2)*tr(rho*A01))
    W[2,1] = real((1/2)*tr(rho*A10))
    W[2,2] = real((1/2)*tr(rho*A11))

    #mana = log2(sum(abs.(W)))
    
    #mana = sum(abs.(W)) - sum(W)
    mana = sum([abs(W[i,j]) for i in 1:2 for j in 1:2]) - sum(W)
    println("flaaaaag:",sum(W)) 

    
   return [W,mana]
    
end



#ba=SpinBasis(1//2)
#idop = dm(spindown(ba)) + dm(spinup(ba))
#Tdm = (1/2)*(idop + (1/3^(1/2))*(sigmax(ba) + sigmay(ba) + sigmaz(ba)))

#rho_q = Tdm

#mr = robustness(rho_q)

#println("test: ",mr)

end
