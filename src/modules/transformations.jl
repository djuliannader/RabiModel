module transformations
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("diagonalization.jl")
include("wigner_eig.jl")
include("dynamics.jl")
using .diagonalization
using .dynamics
using .wigner_eig



function rhodetcat(Hint,rho0,tmax,n,acc)
    times = (0.0,tmax)
    tlist=[0,tmax]
    f(u,p,t) = -im*(Hint)*u + im*u*(Hint) 
    prob = ODEProblem(f,rho0,times)
    sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=0.1,saveat = tlist,
                save_everystep = false)
    rhot = sol.u[length(sol.t)]
    return rhot
end




end
