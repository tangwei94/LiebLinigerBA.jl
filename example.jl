using LiebLinigerBA

# example 1
c = 1.0 # interaction strength
L = 10.0 # system length
N = 10 # particle number

# the ground state in N=10 sector  
ψ = ground_state(c, L, N)
# the chemical potential that makes the ground state stay in the N=10 sector 
μ = get_mu(c, L, N) 
# total ground state energy (with chemical potential μ)
Egs = energy(ψ, μ)

@show c, L, N, μ, Egs

# example 2
c = 1.0 
L = 10.0 
μ = 2.0 
Nmax = 20

# go through different particle numbers
ψs = [ground_state(c, L, Nx) for Nx in 2:Nmax]
Es = [energy(ψi, μ) for ψi in ψs]

# find the ground state
Egs, gs_index = findmin(Es)
ψgs = ψs[gs_index]
Ngs = particle_number(ψgs)

if particle_number(ψgs) == Nmax
    @warn("increase Nmax")
end

@show c, L, μ, Ngs, Egs