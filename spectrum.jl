using LinearAlgebra
using Plots 
using JLD2

using Revise 
using LiebLinigerBA

c = 2.0
L = 20  
N = L

ψ0 = ground_state(c, L, N);

ψps = LLBAState[];
ψms = LLBAState[];
for n in 1:5
    push!(ψps, ground_state(c, L, N+n))
    push!(ψms, ground_state(c, L, N-n))
end

μ = (energy(ψps[1], 0) - energy(ψms[1], 0)) / 2
@show μ
Egs = energy(ψ0, μ)

ph_actions = [[0], [-1], [-2], [-1, -1]];
states = LLBAState[];
ψph1 = ph_excitation(ψ0, [0], [-1]);
ΔE = energy(ψph1, μ) - Egs

for qL in ph_actions
    for qR in ph_actions
        for ψx in [ψ0; ψps; ψms]
            push!(states, ph_excitation(ψx, qL, qR))
        end
    end
end

momenta = momentum.(states)
energies = energy.(states, μ)
scale_E(x::Float64) = (x - Egs) / ΔE
scale_E(ψ::LLBAState, μ::Float64) = (energy(ψ, μ) - Egs) / ΔE

ϕ11 = ph_excitation(ψ0, [0], [-1, -1]);
scale_E(ϕ11, μ)
ϕ2 = ph_excitation(ψ0, [0], [-2]);
scale_E(ϕ2, μ)

ϕp1_11 = ph_excitation(ψps[1], [0], [-1, -1]);
scale_E(ϕp1_11, μ)

ϕm1_11 = ph_excitation(ψms[1], [0], [-1, -1]);
scale_E(ϕm1_11, μ)

plot(momenta, scale_E.(energies), seriestype = :scatter, ylims = [-0.1, 2.5], legend = false)
scatter!([momentum(ϕ11)], [scale_E(ϕ11, μ)], marker = (:xcross, 6, 0.6, :red))
scatter!([momentum(ϕ2)], [scale_E(ϕ2, μ)], marker = (:xcross, 6, 0.6, :blue))

#msk = (abs.(scale_E.(energies) .- 2.0) .< 0.05) .* (momentum.(states) .≈ 4*pi/L)
#ϕ_special = states[msk]
#ϕ_special[1].ns
#scale_E(ϕ_special[1], μ)

# plot arrows to show ph excitations
pnumber(x::LLBAState) = length(x.ns)
colors = cgrad(:roma, 10, categorical = true, scale = :exp) 

N_plt = N
scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)

msk = (pnumber.(states) .== N_plt)
scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)

N_plt = N - 1
scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)

msk = (pnumber.(states) .== N_plt)
scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)

N_plt = N + 1
scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)

msk = (pnumber.(states) .== N_plt)
scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)

jldsave("lieb_liniger_BA_c$(c)_N$(N)_L$(L).jld2"; momenta, energies, Egs)