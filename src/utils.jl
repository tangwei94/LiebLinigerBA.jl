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

function twobody_wave_function(ψ::LLBAState{<:AbstractFloat})
    λ1, λ2 = ψ.quasimomenta
    N = length(ψ.ns)
    L = ψ.L
    c = ψ.c

    if N != 2 
        throw(ArgumentError("two body wavefunction"))
    end

    a, b = exp(im * theta(λ1, λ2, c) / 2), -exp(im * theta(λ2, λ1, c) / 2)    
    @show b/a, exp(-im * λ2 * L) 

    function f(x1::Real, x2::Real)
        if x1 < x2
            return a*exp(im*λ1*x1 + im*λ2*x2) + b*exp(im*λ2*x1 + im*λ1*x2)
        else
            return f(x2, x1)
            #return b*exp(im*λ1*x1 + im*λ2*x2) + a*exp(im*λ2*x1 + im*λ1*x2)
        end
    end

    return f
end