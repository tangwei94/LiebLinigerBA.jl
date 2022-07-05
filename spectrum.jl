using LinearAlgebra
using Plots 
using JLD2

using Revise 
using LiebLinigerBA

c = 8.0
L = 10  
N = 2*L

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
v = ΔE / (2*pi / L)

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

ϕ1 = ph_excitation(ψ0, [0], [-1]);
scale_E(ϕ1, μ)

ϕ11 = ph_excitation(ψ0, [0], [-1, -1]);
scale_E(ϕ11, μ)
ϕ2 = ph_excitation(ψ0, [0], [-2]);
scale_E(ϕ2, μ)

ϕp1_11 = ph_excitation(ψps[1], [0], [-1, -1]);
scale_E(ϕp1_11, μ)

ϕm1_11 = ph_excitation(ψms[1], [0], [-1, -1]);
scale_E(ϕm1_11, μ)

holomorphic_ovlps_gs = [kacmoody(ψi, ψ0, :holomorphic, v) for ψi in states]
anti_holomorphic_ovlps_gs = [kacmoody(ψi, ψ0, :antiholomorphic, v) for ψi in states]
holomorphic_ovlps_m1 = [kacmoody(ψi, ψms[1], :holomorphic, v) for ψi in states]
anti_holomorphic_ovlps_m1 = [kacmoody(ψi, ψms[1], :antiholomorphic, v) for ψi in states]
holomorphic_ovlps_ph1 = [kacmoody(ψi, ϕ1, :holomorphic, v) for ψi in states]
anti_holomorphic_ovlps_ph1 = [kacmoody(ψi, ϕ1, :antiholomorphic, v) for ψi in states]

function tmp_round_off(x)
    if x > N+2
        return N+3
    elseif x<N-2
        return N+2
    else
        return x
    end
end

pnumbers = [tmp_round_off(length(ψi.ns)) for ψi in states]
maximum(pnumbers)
minimum(pnumbers)
23-18

#colors = cgrad(:Blues_3, 100, categorical = true)  
#scatter(momenta, scale_E.(energies), markersize=4, markerstrokewidth=0, markeralpha=0.6, color=colors[end], ylims=(0, 3.5), legend=false)
#plot!(size=(300, 300))
#savefig("ba_spect_full.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(holomorphic_ovlps_gs), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Rgs.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(anti_holomorphic_ovlps_gs), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Lgs.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(holomorphic_ovlps_m1), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Rnp1.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(anti_holomorphic_ovlps_m1), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Lnp1.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(holomorphic_ovlps_ph1), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Re1.pdf")

scatter(momenta, scale_E.(energies), marker_z = norm.(anti_holomorphic_ovlps_ph1), markersize=4, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors[end], markeralpha=0.6, color=colors, ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_Le1.pdf")

colors1 = palette(:tableau_10)
shapes1 = [:diamond, :circle, :star5, :circle, :diamond, :xcross]
msk_gs = [(length(ψi.ns) == N) for ψi in states]
msk_gsp1 = [(length(ψi.ns) == N+1) for ψi in states]
msk_gsp2 = [(length(ψi.ns) == N+2) for ψi in states]
msk_gsp3 = [(length(ψi.ns) == N+3) for ψi in states]
msk_gsm1 = [(length(ψi.ns) == N-1) for ψi in states]
msk_gsm2 = [(length(ψi.ns) == N-2) for ψi in states]
msk_gsm3 = [(length(ψi.ns) == N-3) for ψi in states]
msk_else = [(abs(length(ψi.ns) - N) > 3) for ψi in states]

scatter(momenta[msk_gs], scale_E.(energies[msk_gs]), markersize=4, markershape=:square, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[3], markeralpha=0.6, color=colors1[3], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsp1], scale_E.(energies[msk_gsp1]), markersize=4, markershape=:diamond, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[2], markeralpha=0.6, color=colors1[2], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsm1], scale_E.(energies[msk_gsm1]), markersize=4, markershape=:diamond, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[4], markeralpha=0.6, color=colors1[4], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsp2], scale_E.(energies[msk_gsp2]), markersize=4, markershape=:hexagon, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[5], markeralpha=0.6, color=colors1[5], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsm2], scale_E.(energies[msk_gsm2]), markersize=4, markershape=:hexagon, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[7], markeralpha=0.6, color=colors1[7], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsp3], scale_E.(energies[msk_gsp3]), markersize=4, markershape=:circle, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[8], markeralpha=0.6, color=colors1[8], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_gsm3], scale_E.(energies[msk_gsm3]), markersize=4, markershape=:circle, markerstrokewidth=0.5, markerstrokealpha=0.3, markerstrokecolor=colors1[9], markeralpha=0.6, color=colors1[9], ylims=(0, 3.5), legend=true)
scatter!(momenta[msk_else], scale_E.(energies[msk_else]), markersize=3, markershape=:xcross, markerstrokewidth=1, markerstrokealpha=0.6, markerstrokecolor=colors1[end], markeralpha=0.6, color=colors1[end], ylims=(0, 3.5), legend=true)
plot!(size=(300, 300))
savefig("ba_spect_full.pdf")

#plot(momenta, scale_E.(energies), seriestype = :scatter, ylims = [-0.1, 2.5], legend = false)

#scatter!([momentum(ϕ11)], [scale_E(ϕ11, μ)], marker = (:xcross, 6, 0.6, :red))
#scatter!([momentum(ϕ2)], [scale_E(ϕ2, μ)], marker = (:xcross, 6, 0.6, :blue))
#
##msk = (abs.(scale_E.(energies) .- 2.0) .< 0.05) .* (momentum.(states) .≈ 4*pi/L)
##ϕ_special = states[msk]
##ϕ_special[1].ns
##scale_E(ϕ_special[1], μ)
## plot arrows to show ph excitations
#pnumber(x::LLBAState) = length(x.ns)
#colors = cgrad(:roma, 10, categorical = true, scale = :exp) 
#
#N_plt = N
#scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)
#
#msk = (pnumber.(states) .== N_plt)
#scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)
#
#N_plt = N - 1
#scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)
#
#msk = (pnumber.(states) .== N_plt)
#scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)
#
#N_plt = N + 1
#scatter(momenta, scale_E.(energies), marker = (:circle, 4, 0.5, colors[end]), ylims = [-0.1, 2.5], legend = false)
#
#msk = (pnumber.(states) .== N_plt)
#scatter!(momenta[msk], scale_E.(energies[msk]), marker= (:circle, 4, 0.9, colors[abs(N_plt - N) + 1]), legend=false)
#
#jldsave("lieb_liniger_BA_c$(c)_N$(N)_L$(L).jld2"; momenta, energies, Egs)