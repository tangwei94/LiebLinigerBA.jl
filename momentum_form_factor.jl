using LinearAlgebra
using Plots

using Revise
using LiebLinigerBA

c = 16
L = 10
N = L
ψ = ground_state(c, L, N)
ϕ0 = ph_excitation(ψ, [-1], [-1])
ϕ1 = ph_excitation(ψ, [0], [-1])
ϕ2 = ph_excitation(ψ, [-1], [-2])

ψis = [intermediate_states(c, L, N-1, Δp) for Δp in 0:12];
@show length.(ψis)

function momentum_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat})
    result_j = [] 
    result = 0
    for ψis_j in ψis
        for ψi in ψis_j
            ln_tmpR, phase_tmpR = ln_ψ0_form_factor(ψi, ψ2)
            ln_tmpL, phase_tmpL = ln_ψ0_form_factor(ψi, ψ1)
            result += (momentum(ψ2) - momentum(ψi)) * exp(ln_tmpL + ln_tmpR) * (phase_tmpR * phase_tmpL')
        end
        @show result
        push!(result_j, result)
    end
    return result_j
end

p10_j = momentum_form_factor(ϕ1, ϕ0)
p20_j = momentum_form_factor(ϕ2, ϕ0)

lnnorm_ρ10, phase_ρ10 = ln_ρ0_form_factor(ϕ1, ϕ0)
ρ10 = exp(lnnorm_ρ10) * phase_ρ10
lnnorm_ρ20, phase_ρ20 = ln_ρ0_form_factor(ϕ2, ϕ0)
ρ20 = exp(lnnorm_ρ20) * phase_ρ20

p10_j
p10_j[end] / ρ10
p20_j
p20_j[end] / ρ20

plot(Vector(0:12), norm.(result_j))