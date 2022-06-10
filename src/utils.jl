"""
    ph_labels(n_particles::Integer, Δp::Integer)

    Generate a collection of labels for all possible particle-hole excitations:
    - `n_particles`: total number of particles in the branch (L/R sector)
                     serves as an upper limit to the number of ph excitations. 
    - `Δp`: the total momentum increase 
"""
function ph_labels(n_particles::Integer, Δp::Integer)
    if Δp == 0
        return [[0]]
    end

    ph_excitations = []
    for num_devide in 1:min(n_particles, Δp)
        tmp = collect(partitions(Δp, num_devide)) .* (-1)
        append!(ph_excitations, tmp)
    end
    return ph_excitations
end

function ph_labels(n_particles::Integer, Δpmin::Integer, Δpmax::Integer)
    ph_excitations = []

    for Δp in Δpmin:Δpmax
        append!(ph_excitations, ph_labels(n_particles, Δp))
    end
    return ph_excitations
end

"""
    Generate a series of intermediate states for given total particle number n.
"""
function intermediate_states(c::Real, L::Real, n::Integer, Δp::Integer)
    ψis = []

    gs_sectors = [(n-ix, ix) for ix in 1:n]
    for gs_label in gs_sectors 
        nl, nr = gs_label 
        ψi_gs = ground_state(c, L, nl, nr)

        for l_ph in ph_labels(nl, Δp)
            for r_ph in ph_labels(nr, 0, Δp)
                ψi = ph_excitation(ψi_gs, l_ph, r_ph)
                push!(ψis, ψi)
            end
        end

        for l_ph in ph_labels(nl, 0, Δp-1)
            for r_ph in ph_labels(nr, Δp)
                ψi = ph_excitation(ψi_gs, l_ph, r_ph)
                push!(ψis, ψi)
            end
        end
    end
    return ψis
end

function get_mu(c::Real, L::Real, N::Integer)
    E1 = energy(ground_state(c, L, N-1), 0)
    E2 = energy(ground_state(c, L, N+1), 0)
    return (E2 - E1) / 2
end

function Luttinger_parameter(c::Real, L::Real, N::Integer)
    μ = get_mu(c, L, N)

    ψ0 = ground_state(c, L, N)
    E0 = energy(ψ0, μ)
    ψ1 = ph_excitation(ψ0, [0], [-1])
    E1 = energy(ψ1, μ)

    ψ_1p = ground_state(c, L, N+1)
    E_1p = energy(ψ_1p, μ)

    return (E1 - E0) / (E_1p - E0) / 4
end

function twobody_wave_function(ψ::LLBAState{<:AbstractFloat})
    λ1, λ2 = ψ.quasimomenta
    N = length(ψ.ns)
    L = ψ.L
    c = ψ.c

    if N != 2 
        throw(ArgumentError("two body wavefunction"))
    end

    a, b = exp(im * theta(λ1, λ2, c) / 2), -exp(im * theta(λ2, λ1, c) / 2)    

    function f(x1::Real, x2::Real)
        y1, y2 = sort([x1, x2])
        return a*exp(im*λ1*y1 + im*λ2*y2) + b*exp(im*λ2*y1 + im*λ1*y2)
    end

    return f
end

function threebody_wave_function(ψ::LLBAState{<:AbstractFloat})

    λs = ψ.quasimomenta
    N, L, c = length(ψ.ns), ψ.L, ψ.c

    if N!=3
        throw(ArgumentError("three body wavefunction"))
    end

    perms = collect(permutations([1,2,3]))
    λs_perms = [λs[perm] for perm in perms]
    parities = levicivita.(perms)

    coeffs = [η * exp(0.5*im * theta(Pλs[1], Pλs[2], c) + 0.5*im * theta(Pλs[1], Pλs[3], c) + 0.5*im * theta(Pλs[2], Pλs[3], c)) for (Pλs, η) in zip(λs_perms, parities)]

    function f(x1::Real, x2::Real, x3::Real)
        ys = sort([x1, x2, x3])
        return sum(coeffs .* [exp(sum(Pλs .* ys .* im)) for Pλs in λs_perms])
    end
    return f
end
