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

θs = 0:0.01*π:2*π
ellipse_xs = 3 .+ 0.4 .* cos.(θs)
ellipse_ys = 3 .+ 0.4 .* sin.(θs)

function plot_spect(ax, ψi, ovlpsR, ovlpsL)
    msk_Rp = real.(ovlpsR) .> 0.6
    msk_Lp = real.(ovlpsL) .> 0.6
    msk_Rm = real.(ovlpsR) .< -0.6
    msk_Lm = real.(ovlpsL) .< -0.6
    msk_i = isapprox.(energies, energy(ψi)) .& isapprox.(momenta .+ 0.1, momentum(ψi)+ 0.1)  
    msk_0 = .~(msk_Rp .| msk_Rm .| msk_Lp .| msk_Lm .| msk_i)

    colorp = :royalblue
    colorm = :mediumpurple

    sc_main = scatter!(ax, momenta[msk_0] .* L ./ (2*pi), scale_E.(energies[msk_0]), color=:darkorange, marker=:circle, markersize=10)
    sc_mappedRp = scatter!(ax, momenta[msk_Rp] .* L ./ (2*pi), scale_E.(energies[msk_Rp]), color=colorp, marker=:diamond, markersize=10, label=L"\text{positive overlap with } J_n |\psi_\mathrm{i}\rangle")
    sc_mappedRm = scatter!(ax, momenta[msk_Rm] .* L ./ (2*pi), scale_E.(energies[msk_Rm]), color=colorm, marker=:diamond, markersize=10, label=L"\text{negative overlap with } J_n |\psi_\mathrm{i}\rangle")
    sc_mappedLp = scatter!(ax, momenta[msk_Lp] .* L ./ (2*pi), scale_E.(energies[msk_Lp]), color=colorp, marker=:rect, markersize=10, label=L"\text{positive overlap with } \bar{J}_n |\psi_\mathrm{i}\rangle")
    sc_mappedLm = scatter!(ax, momenta[msk_Lm] .* L ./ (2*pi), scale_E.(energies[msk_Lm]), color=colorm, marker=:rect, markersize=10, label=L"\text{negative overlap with } \bar{J}_n |\psi_\mathrm{i}\rangle")
    return sc_main, sc_mappedRp, sc_mappedRm, sc_mappedLp, sc_mappedLm
end

#font1 = Makie.to_font("/home/wtang/.local/share/fonts/STIXTwoText-Regular.otf")
fig = Figure(backgroundcolor = :white, fontsize=18, resolution= (600, 425))#, fonts=(; regular=font1))
gf = fig[1, 1] = GridLayout() 
gl = fig[2, 1] = GridLayout()

ax1 = Axis(gf[1, 1], 
        xlabel = L"p L / 2\pi",
        ylabel = L"\text{rescaled } \Delta E", 
        xticks = -4:1:4,
        yticks = 0:1:4,
        )

ylims!(ax1, (-0.15, 4.5))

sc1_orig = scatter!(ax1, [0], [0], color=:gray50, marker=:star5, markersize=15)
sc1_main, sc1_mappedRp, sc1_mappedRm, sc1_mappedLp, sc1_mappedLm = plot_spect(ax1, ψ0, holomorphic_ovlps_gs, anti_holomorphic_ovlps_gs)
lines!(ax1, ellipse_xs, ellipse_ys, color=:gray66)

@show fig

ax2 = Axis(gf[1, 2], 
        xlabel = L"p L / 2\pi",
        #ylabel = L"\text{rescaled } \Delta E", 
        xticks = -4:1:4,
        yticks = 0:1:4,
        )

ylims!(ax2, (-0.15, 4.5))

sc2_main, sc2_mappedRp, sc2_mappedRm, sc2_mappedLp, sc2_mappedLm  = plot_spect(ax2, ψ1, holomorphic_ovlps_1, anti_holomorphic_ovlps_1l)
sc2_orig = scatter!(ax2, [1], [scale_E(energy(ψ1))], color=:gray50, marker=:star5, markersize=15)
lines!(ax2, ellipse_xs, ellipse_ys, color=:gray66)

@show fig

for (label, layout) in zip(["(a)", "(b)"], [gf[1, 1], gf[1, 2]])
    Label(layout[1, 1, TopLeft()], label, 
    padding = (0, 5, -5, 0), 
    halign = :right
    )
end

Label(gf[1,1][1, 1, TopLeft()], L"|\psi_\mathrm{i}\rangle",
    padding = (0, -310, -470, 0))
Label(gf[1,2][1, 1, TopLeft()], L"|\psi_\mathrm{i}\rangle",
    padding = (0, -365, -370, 0))

#axislegend(ax1, position=:lb, framevisible=false)

leg = Legend(gl[1,1], ax1, orientation=:horizontal, framecolor=:lightgrey, labelsize=16)
leg.nbanks = 2

@show fig
Makie.save("fig-ba-spect.pdf", fig)
