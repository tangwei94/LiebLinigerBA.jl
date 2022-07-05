using LinearAlgebra
using Plots

using Revise
using LiebLinigerBA

c = 8
L = 10
N = 20
γ = c * L / N
K = Luttinger_parameter(c, L, N)
v = velocity(c, L, N)
@show c, L, N, γ, v, K

ψ = ground_state(c, L, N);
ϕ0 = ph_excitation(ψ, [-2], [0]);
ϕ1 = ph_excitation(ψ, [0], [0]);
ϕ2 = ph_excitation(ψ, [-2], [-1, -1]);

#ψis = [intermediate_states(c, L, N-1, Δp) for Δp in 0:12];
#@show length.(ψis)

#function momentum_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat})
#    result_j = [] 
#    result = 0
#    for ψis_j in ψis
#        for ψi in ψis_j
#            ln_tmpR, phase_tmpR = ln_ψ0_form_factor(ψi, ψ2)
#            ln_tmpL, phase_tmpL = ln_ψ0_form_factor(ψi, ψ1)
#            result += (momentum(ψ2) + momentum(ψ1) - 2*momentum(ψi)) * exp(ln_tmpL + ln_tmpR) * (phase_tmpR * phase_tmpL')
#        end
#        #@show result
#        push!(result_j, result)
#    end
#    return result_j
#end

#j10_j = momentum_form_factor(ϕ1, ϕ0) 
lnnorm_j10, phase_j10 = ln_ρ0_form_factor(ϕ1, ϕ0; target=:current)
j10 = exp(lnnorm_j10) * phase_j10
lnnorm_j20, phase_j20 = ln_ρ0_form_factor(ϕ2, ϕ0; target=:current)
j20 = exp(lnnorm_j20) * phase_j20

lnnorm_ρ10, phase_ρ10 = ln_ρ0_form_factor(ϕ1, ϕ0)
ρ10 = exp(lnnorm_ρ10) * phase_ρ10
lnnorm_ρ20, phase_ρ20 = ln_ρ0_form_factor(ϕ2, ϕ0)
ρ20 = exp(lnnorm_ρ20) * phase_ρ20

@show ρ10 / j10 / (L/N/pi/2) / K
@show ρ20 / j20 / (L/N/pi/2) / K
@show ρ10 / j10 * v 
@show ρ20 / j20 * v 
@show 1 / L * N * pi *2 / K  - v

@show ρ10 * (2*pi*N/L) / K + j10
@show ρ10 * (2*pi*N/L) / K - j10
@show ρ20 * (2*pi*N/L) / K + j20
@show ρ20 * (2*pi*N/L) / K - j20

