module QPT
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
import troterization
import wigner_eig
import stat
export amplitud


function QPTexpectation(Nmax::Int64,nn::Int64,r,om,lambda,delta,chi,nu,eta,psi)
         etar=0
	 etaint=0.2
	 expval=[]
	 vdiag = [(i) + 0*im for i in 0:Nmax for j in -1:0]
	 numop=Array(Diagonal(vdiag))
         for i in 1:60
	   floquet=troterization.troter(Nmax,nn,r,om,lambda,delta,chi,nu,etar,psi)
	   vecordered= stat.orderinfvec(Nmax,om,r,lambda,delta,nn,nu,chi,etar,psi)
	   #ham=diagonalization.hamiltonian(Nmax,om,r,lambda,delta,etar,psi)
	   evfloquet = eigvecs(floquet)
	   #evfloquet = eigvecs(ham)
	   ket = [evfloquet[i,vecordered[1]] for i in 1:2(Nmax+1)]
	   #ket = [evfloquet[i,1] for i in 1:2(Nmax+1)]
	   bra=Array{Complex{Float64}}(undef,1,length(ket)) 
           for k in 1:length(ket)
              bra[1,k]=conj(ket[k])
           end
	   expvalinst= bra * numop * ket
           append!(expval,expvalinst[1])
	 etar=etar + etaint
         end
	 return expval
end

function QPTexpectation_rmp(Nmax::Int64,nn::Int64,r,om,omega,eta,psi,gamma)
         etar=0.0
	 etaint=0.5
	 expval=[]
	 vdiag = [(i) + 0*im for i in 0:Nmax for j in -1:0]
	 numop=Array(Diagonal(vdiag))
         for i in 1:40
	   floquet=troterization.troter2(Nmax,nn,r,om,gamma,omega,etar,psi)
	   vecordered= stat.orderinfvec2(Nmax,om,r,gamma,omega,etar,psi,nn)
	   evfloquet = eigvecs(floquet)
	   #evfloquet = eigvecs(ham)
	   ket = [evfloquet[i,vecordered[1]] for i in 1:2(Nmax+1)]
	   bra=Array{Complex{Float64}}(undef,1,length(ket)) 
           for k in 1:length(ket)
              bra[1,k]=conj(ket[k])
           end
	   expvalinst= bra * numop * ket
           append!(expval,expvalinst[1])
	 etar=etar + etaint
	 println(i)
         end
	 println(etar)
	 #msj=wigner_eig.wigner_driven2(Nmax,om,r,gamma,omega,etar,psi,500,1,10)
	 return expval
end

n=50
nn=500
om=1.0
r=50.0
lambda=1.0
delta=0.0
eta=5.0
psi=0.5
omega=1.0
gamma=1.0
b=0.01
nu=1.0



#QPTlist = QPTexpectation(n,nn,r,om,lambda,delta,b,nu,eta,psi)
QPTlist = QPTexpectation_rmp(n,nn,r,om,omega,eta,psi,gamma)
println(QPTlist)


end
