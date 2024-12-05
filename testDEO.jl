using DifferentialEquations

times=(0.0,1.0)

H=[0.0 0.0 0.0 0.0; 0.0 0.0 1.0*im 0.0; 0.0 1.0*im 0.0 0.0; 0.0 0.0 0.0 0.0]
u0 = [0.0+0.0*im 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im; 0.0+0.0*im 1.0+0.0*im 0.0+0.0*im 0.0+0.0*im; 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im; 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im]
f(u,p,t) = -1.0*(H*u-u*H)
prob = ODEProblem(f,u0,times)
sol = solve(prob)

println(sol(0.5))
println(H)


