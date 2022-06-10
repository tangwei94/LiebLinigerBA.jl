using LinearAlgebra

using Revise
using LiebLinigerBA

c = 32
L = 10
N = L
ψ = ground_state(c, L, N)
ϕ = ph_excitation(ψ, [-2], [-1])

ψis = [intermediate_states(c, L, N-1, Δp) for Δp in 0:12];
@show length.(ψis)

result_j = [] 
result = 0
for ψis_j in ψis
    for ψi in ψis_j
        ln_tmpR, phase_tmpR = ψ0_form_factor(ψi, ψ)
        ln_tmpL, phase_tmpL = ψ0_form_factor(ψi, ϕ)
        result += (momentum(ψ) - momentum(ψi)) * exp(ln_tmpL + ln_tmpR) * (phase_tmpR * phase_tmpL')
    end
    @show result
    push!(result_j, result)
end