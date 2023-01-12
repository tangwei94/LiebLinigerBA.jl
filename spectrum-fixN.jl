using LinearAlgebra
using LaTeXStrings
using CairoMakie

using Revise 
using LiebLinigerBA

c = 1.
L = 64.
N = Int(L) 
ρ0 = N / L 

# the ground state
ψ0 = ground_state(c, L, N);
ψ1 = ph_excitation(ψ0, [0], [-1]);
ψ1l = ph_excitation(ψ0, [-1], [0]);
Egs = energy(ψ0);
v, K = v_and_K(c, L, N);

ΔE = energy(ψ1) - Egs
scale_E(x::Float64) = (x - Egs) / ΔE
scale_E(ψ::LLBAState) = (energy(ψ) - Egs) / ΔE

# excited states
ph_actions = [[0], [-1], [-2], [-1, -1], [-3], [-2, -1], [-1, -1, -1], [-4], [-3, -1], [-2, -2], [-2, -1, -1], [-1, -1, -1, -1]];
states = LLBAState[];
for qL in ph_actions
    for qR in ph_actions
        push!(states, ph_excitation(ψ0, qL, qR))
    end
end

msk = scale_E.(energy.(states)) .< 4.5
states = states[msk];

momenta = momentum.(states)
energies = energy.(states)

holomorphic_ovlps_gs = [kacmoody(ψi, ψ0, :holomorphic, v, K) for ψi in states]
anti_holomorphic_ovlps_gs = [kacmoody(ψi, ψ0, :antiholomorphic, v, K) for ψi in states]
holomorphic_ovlps_1 = [kacmoody(ψi, ψ1, :holomorphic, v, K) for ψi in states]
anti_holomorphic_ovlps_1l = [kacmoody(ψi, ψ1, :antiholomorphic, v, K) for ψi in states]

function plot_spect(ax, ψi, ovlpsR, ovlpsL)
    msk_R = norm.(ovlpsR) .> 0.6
    msk_L = norm.(ovlpsL) .> 0.6
    msk_i = isapprox.(energies, energy(ψi)) .& isapprox.(momenta .+ 0.1, momentum(ψi)+ 0.1)  
    msk_0 = .~(msk_R .| msk_L .| msk_i)
    
    colorR = :royalblue
    colorL = :mediumpurple

    sc_main = scatter!(ax, momenta[msk_0] .* L ./ (2*pi), scale_E.(energies[msk_0]), color=:darkorange, marker=:circle, markersize=10)
    sc_mappedR = scatter!(ax, momenta[msk_R] .* L ./ (2*pi), scale_E.(energies[msk_R]), color=colorR, marker=:diamond, markersize=10, label=L"\text{mapped from } |\psi_i\rangle \text{ with } J_{n}")
    sc_mappedL = scatter!(ax, momenta[msk_L] .* L ./ (2*pi), scale_E.(energies[msk_L]), color=colorL, marker=:diamond, markersize=10, label=L"\text{mapped from } |\psi_i\rangle \text{ with } \bar{J}_{n}")
    return sc_main, sc_mappedR, sc_mappedL
end

fig = Figure(backgroundcolor = :white, fontsize=18, resolution= (600, 400))
gf = fig[1, 1] = GridLayout() 
gl = fig[2, 1] = GridLayout()

ax1 = Axis(gf[1, 1], 
        xlabel = L"p L / 2\pi",
        ylabel = L"\text{rescaled } \Delta E", 
        xticks = -4:1:4,
        yticks = 0:1:4,
        )

ylims!(ax1, (-0.15, 4.5))

sc1_orig = scatter!(ax1, [0], [0], color=:gray50, marker=:star5, markersize=15, label=L"|\psi_{i}\rangle")
sc1_main, sc1_mappedR, sc1_mappedL = plot_spect(ax1, ψ0, holomorphic_ovlps_gs, anti_holomorphic_ovlps_gs)

@show fig

ax2 = Axis(gf[1, 2], 
        xlabel = L"p L / 2\pi",
        #ylabel = L"\text{rescaled } \Delta E", 
        xticks = -4:1:4,
        yticks = 0:1:4,
        )

ylims!(ax2, (-0.15, 4.5))

sc2_main, sc2_mappedR, sc2_mappedL = plot_spect(ax2, ψ1, holomorphic_ovlps_1, anti_holomorphic_ovlps_1l)
sc2_orig = scatter!(ax2, [1], [scale_E(energy(ψ1))], color=:gray50, marker=:star5, markersize=15)

@show fig

for (label, layout) in zip(["(a)", "(b)"], [gf[1, 1], gf[1, 2]])
    Label(layout[1, 1, TopLeft()], label, 
    padding = (0, 5, -5, 0), 
    halign = :right
    )
end

#axislegend(ax1, position=:lb, framevisible=false)

leg = Legend(gl[1,1], ax1, orientation=:horizontal, framecolor=:lightgrey, labelsize=16)

@show fig
Makie.save("fig-ba-spect.pdf", fig)
