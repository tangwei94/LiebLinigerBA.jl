using LiebLinigerBA

c = 1.0 # interaction strength
L = 10.0 # system length
N = 10 # particle number

# the ground state in N=10 sector  
ψ = ground_state(c, L, N)
# the chemical potential that makes the ground state stay in the N=10 sector 
μ = get_mu(c, L, N) 
# ground state energy (with chemical potential μ)
Egs = energy(ψ, μ)

@show c, L, N, μ, Egs