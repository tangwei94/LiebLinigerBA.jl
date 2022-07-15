"""
    struct LLBAState{T<:AbstractFloat}

    represent a Lieb-Liniger Bethe ansatz state.
"""
struct LLBAState{T<:AbstractFloat}
    ns::Vector{<:Rational}
    quasimomenta::Vector{T}
    c::Real
    L::Real

    function LLBAState(ns::Vector{<:Rational}, quasimomenta::Vector{T}, c::Real, L::Real) where T<:AbstractFloat
        # sort input ns and quasimomenta
        ns = sort(ns)
        quasimomenta = sort(quasimomenta)

        if all(denominator.(ns) .== 1) 
            length(ns) % 2 == 1 || throw(ArgumentError("integer ns: particle number should be odd")) 
        elseif all(denominator.(ns) .== 2) 
            length(ns) % 2 == 0 || throw(ArgumentError("half-integer ns: particle number should be even")) 
        else
            throw(ArgumentError("ns can only be all-half-integer or all-integer"))
        end

        if length(ns) != length(quasimomenta)
            throw(ArgumentError("ns should have equal length with quasimomenta"))
        end

        if prod(ns[2:end] - ns[1:end-1]) ≈ 0 
            throw(ArgumentError("input ns cannot contain equal numbers"))
        end

        # verify the BA equation
        for (λ, n) in zip(quasimomenta, ns)
            if abs(λ * L - 2*pi*n - sum(theta.(λ, quasimomenta, c)) ) > 1e-8
                throw(ArgumentError("input ns and quasimomenta don't satisfy the BA equations"))
            end
        end

        return new{T}(ns, quasimomenta, c, L)
    end
end

function LLBAState(ns::Vector{<:Rational}, c::Real, L::Real)
    quasimomenta = Float64.(copy(ns))
    BA_solve!(ns, quasimomenta, c, L)
    return LLBAState(ns, quasimomenta, c, L)
end

LLBAState(ns::Vector{<:Integer}, c::Real, L::Real) = LLBAState(Rational.(ns), c, L)
LLBAState(ns::Vector{<:Integer}, quasimomenta::Vector{<:AbstractFloat}, c::Real, L::Real) = LLBAState(Rational.(ns), quasimomenta, c, L)

energy(ψ::LLBAState{<:AbstractFloat}, μ::Real) = sum(ψ.quasimomenta .^ 2) - μ * length(ψ.ns)
particle_number(ψ::LLBAState{<:AbstractFloat}) = length(ψ.ns)
momentum(ψ::LLBAState{<:AbstractFloat}) = sum(ψ.quasimomenta) 

"""
    theta(λ1::Real, λ2::Real, c::Real)

    calculate θ(λ1, λ2). `λ1` and `λ2` are quasimomenta; `c` is the parameter in the Hamiltonian. 
"""
theta(λ1::Real, λ2::Real, c::Real) = -2 * atan((λ1 - λ2) / c)

"""
    gtheta(λ1::Real, λ2::Real, c::Real)

    calculate θ'(λ1, λ2). `λ1` and `λ2` are quasimomenta; `c` is the parameter in the Hamiltonian. 
"""
gtheta(λ1::Real, λ2::Real, c::Real) = -2 * c / ((λ1 - λ2)^2 + c^2)

"""
    BA_solve!(ns::Vector{<:Rational}, quasimomenta::Vector{<:AbstractFloat}, c::Real, L::Real)

    Given `ns`, solve the BA equations by minimizing the Yang-Yang action. 
    The vector `quasimomenta` will be modified. The input `quasimomenta` will be used as an initialization. 
"""
function BA_solve!(ns::Vector{<:Rational}, quasimomenta::Vector{<:AbstractFloat}, c::Real, L::Real)
    N = length(ns)
    function YangYang_action_gh(λs)
        grad = L .* λs - 2 * pi .* ns - vec(sum(theta.(λs', λs, c), dims=1))
        hess = gtheta.(λs', λs, c)
        hess += Diagonal(L .* ones(N) - vec(sum(hess, dims=1)))  
        return grad, hess
    end

    gnorm = 1
    while gnorm >= 1e-8
        grad, hess = YangYang_action_gh(quasimomenta)
        gnorm = norm(grad)
        quasimomenta .+= cholesky(hess) \ (-grad)
    end
    return quasimomenta
end

"""
    update_state(ψ::LLBAState{T}, ns::Vector{<:Rational}) where {T}
    update_state(ψ::LLBAState{T}, ns::Vector{<:INteger}) where {T}

    update a state `ψ` with new `ns`. 
    The quasimomenta of `ψ` will be used as initialization.
    Return the updated new state. the original `ψ` won't be changed.
"""
function update_state(ψ::LLBAState{T}, ns::Vector{<:Rational}) where T
    quasimomenta = copy(ψ.quasimomenta)
    BA_solve!(ns, quasimomenta, ψ.c, ψ.L)
    return LLBAState(ns, quasimomenta, ψ.c, ψ.L) 
end
update_state(ψ::LLBAState{T}, ns::Vector{<:Integer}) where T = update_state(ψ, Rational.(ns))

"""
    ground_state(c::Real, L::Real, N::Integer)

    Obtain the ground state in the sector of particle number `N`.
"""
function ground_state(c::Real, L::Real, N::Integer)
    ns = Vector((-N+1:2:N) .// 2)
    return LLBAState(ns, c, L)
end

"""
    ground_state(c::Real, L::Real, nL::Integer, nR::Integer)

    Obtain the ground state in the following sector: 
    there are `nL` particles in the left branch, and `nR` particles in the right branch.
    If the total particle number is odd, then `nR` also takes account of the quantum number at zero, i.e. nR (in program) = nR (actual) + 1. 
"""
function ground_state(c::Real, L::Real, nL::Integer, nR::Integer)
    if (nL + nR) % 2 == 1
        ns = -nL:nR-1
    else
        ns = (-nL+1:nR) .- 1//2
    end
    
    return LLBAState(Vector(ns), c, L)
end

"""
    ph_excitation(ψ::LLBAState{<:AbstractFloat}, qLs::Vector{Integer}, qRs::Vector{Integer}) 

    Generate particle-hole excitations on a given state `ψ`. 
    `qLs` and `qRs` represent the particle-hole excitations on the left and right branch.
    A integer number `q` in `qLs`/`qRs` corresponds to `j_{q}` operator in U1 Kac-Moody algebra: q<0 -> generate p-h excitation; q>0 -> annihilate p-h excitation.  
"""
function ph_excitation(ψ::LLBAState{<:AbstractFloat}, qLs::Vector{<:Integer}, qRs::Vector{<:Integer}) 
    ns = copy(ψ.ns)
    for (ix, qL) in enumerate(qLs) 
        ns[ix] -= -qL
    end
    for (ix, qR) in enumerate(qRs)
        ns[end-(ix-1)] += -qR
    end
    
    if ns != sort(ns)
        throw(ArgumentError("the ph-excitation should not reorder ns"))
    end

    if min(length(qLs), length(qRs)) > length(ns) / 2
        @warn "number of particles too small compared to number of p-h exciations"
    end 

    return update_state(ψ, ns)
end

