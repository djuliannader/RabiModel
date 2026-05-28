module tt
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/diagonalization.jl")
include("modules/reading.jl")
include("modules/stat.jl")
include("modules/troterization.jl")
#include("modules/dynamics.jl")
include("modules/wigner_eig.jl")
using .diagonalization
using .reading
using .stat
using .troterization
#using .dynamics
using .wigner_eig


j=1
i=10

println(i*j)




end
