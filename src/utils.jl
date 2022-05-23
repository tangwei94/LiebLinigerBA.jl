#gs_ns_cfg(N::integer) = Vector((-N+1:2:N) .// 2)

function ground_state(c::Real, L::Real, N::Integer)
    ns = Vector((-N+1:2:N) .// 2)
    return LLBAState(ns, c, L)
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

