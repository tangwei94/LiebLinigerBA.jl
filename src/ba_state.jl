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

"""
    ln_gaudin_norm(ψ::LLBAState{<:AbstractFloat})

    Compute the log of the Gaudin norm of a LLBAState `log[sqrt(<ψ|ψ>)]`. 
"""
function ln_gaudin_norm(ψ::LLBAState{<:AbstractFloat})
    λs = ψ.quasimomenta
    c, L = ψ.c, ψ.L
    N = length(ψ.ns)
    Yang_hess = gtheta.(λs', λs, c)
    Yang_hess += Diagonal(L .* ones(N) - vec(sum(Yang_hess, dims=1)))  
    result = log(det(Yang_hess)) + N * log(c)

    for k in 1:N
        for j in (k+1):N
            λjk = λs[j] - λs[k]
            result += log((λjk^2 + c^2) / λjk^2)
        end
    end
    return 0.5*result
end

"""
    ln_ρ0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0)

    form_factor of the density operator <ψ1|ψ†(0)ψ(0)|ψ2>. 
    ψ1 and ψ2 will be normalized.
    return the log of the norm as well as the total phase.
    Ref: J. De Nardis, M. Panfil, J. Stat. Mech. 2015, P02019 (2015)

    parameter `p` is a free parameter in the formula which will not affect the result.
"""
function ln_ρ0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0)

    if ψ1.ns == ψ2.ns
        return log(length(ψ1.ns) / ψ1.L) 
    end

    λs = ψ1.quasimomenta
    μs = ψ2.quasimomenta
    c, N = ψ1.c, length(ψ1.ns)

    Vps = ones(ComplexF64, N)
    Vms = ones(ComplexF64, N)
    V0s = ones(Float64, N)
    for j in 1:N
        for (λ, μ) in zip(λs, μs)
            Vps[j] *= (μ - λs[j] + im*c) / (λ - λs[j] + im*c)
            Vms[j] *= (μ - λs[j] - im*c) / (λ - λs[j] - im*c)
            if (abs(λ - λs[j]) > 1e-8)
                V0s[j] *= (μ - λs[j]) / (λ - λs[j])
            end
        end
    end

    Vp_0, Vm_0 = 1, 1
    for (λ, μ) in zip(λs, μs)
        Vp_0 *= (μ - p + im*c) / (λ - p + im*c)
        Vm_0 *= (μ - p - im*c) / (λ - p - im*c)
    end

    Id_p_U = Matrix{ComplexF64}(I, (N, N))
    for j in 1:N
        factor = im * (μs[j] - λs[j]) / (Vps[j] - Vms[j]) * V0s[j]
        for k in 1:N
            Id_p_U[j, k] += factor * (gtheta(p, λs[k], c) - gtheta(λs[j], λs[k], c))
        end
    end

    result = det(Id_p_U) / (Vp_0 - Vm_0) * sum(μs .- λs) * prod(Vps .- Vms)
    result_ln_norm = log(norm(result))
    result_angle = angle(result)

    for j in 1:N
        for k in 1:N 
            term_jk = (λs[j] - λs[k] + im*c) / (μs[j] - λs[k])
            result_ln_norm += log(norm(term_jk))
            result_angle += angle(result)
        end
    end

    ln_norm1, ln_norm2 = ln_gaudin_norm(ψ1), ln_gaudin_norm(ψ2)
    result_ln_norm -= ln_norm1 + ln_norm2
    phase = exp(im * result_angle)

    return result_ln_norm, phase
end

"""
    form_factor of the local interaction operator <ψ1|ψ†(0)ψ†(0)ψ(0)ψ(0)|ψ2>. 
    ψ1 and ψ2 will be normalized.
    
    L. Piroli, P. Calabrese, J. Phys. A: Math. Theor. 48, 454002 (2015)
"""
function ln_I0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0, s::Real=0)
    c = ψ1.c
    N = length(ψ1.ns)

    μs = ψ1.quasimomenta
    λs = ψ2.quasimomenta

    function fV(λ, vu, vd)
        result = 1
        for ix in 1:N 
            result *= (vu[ix] - λ) / (vd[ix] - λ)
        end
        return result
    end
    fpV(λ) = fV(λ, μs .+ im*c, λs .+ im*c)
    fmV(λ) = fV(λ, μs .- im*c, λs .- im*c)

    Vps = [fpV(λ) for λ in λs]
    Vms = [fmV(λ) for λ in λs]
    V0s = [fV(λ, λs .- im*c, μs) for λ in λs]

    Id_p_U = Matrix{ComplexF64}(I, (N, N))
    for j in 1:N
        factor = im / (Vps[j] - Vms[j]) * fV(λs[j], μs, [λs[1:j-1]; λs[j] + 1; λs[j+1:end]])
        for k in 1:N
            Id_p_U[j, k] += factor * (-gtheta(λs[j], λs[k], c) - gtheta(p, λs[k], c) * gtheta(s, λs[j], c))
        end
    end

    Pλ, Pμ = sum(λs), sum(μs)
    Eλ, Eμ = sum(λs .^ 2), sum(μs .^ 2)
    Qλ, Qμ = sum(λs .^ 3), sum(μs .^ 3)
    J = (Pλ - Pμ)^4 - 4 * (Pλ - Pμ) * (Qλ - Qμ) + 3 * (Eλ - Eμ)^2

    result = det(Id_p_U) / (fpV(p) - fmV(p)) / (fpV(s) - fmV(s))
    result_ln_norm = log(norm(result))
    result_angle = angle(result)

    for j in 1:N
        vj = (Vps[j] - Vms[j]) * V0s[j]
        result_ln_norm += log(norm(vj))
        result_angle += angle(vj)
    end

    result_ln_norm += log(norm(J/6/c)) - ln_gaudin_norm(ψ1) - ln_gaudin_norm(ψ2)

    return result_ln_norm, exp(im*result_angle) * (-1)^N * sign(J) 
end

"""
    form_factor of the field operator <ψ1|ψ(0)|ψ2>. ψ1 and ψ2 will be normalized.

    L. Piroli, P. Calabrese, J. Phys. A: Math. Theor. 48, 454002 (2015)
"""
function ψ0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0, s::Real=0)
    c = ψ1.c
    N = length(ψ1.ns)

    μs = ψ1.quasimomenta
    λs = ψ2.quasimomenta

    function fV(λ, vu, vd)
        result = 1
        for ix in 1:N 
            result *= (vu[ix] - λ) / (vd[ix] - λ)
        end
        result /= (vd[N+1] - λ)
        return result
    end
    fpV(λ) = fV(λ, μs .+ im*c, λs .+ im*c)
    fmV(λ) = fV(λ, μs .- im*c, λs .- im*c)

    Vps = [fpV(λ) for λ in λs]
    Vms = [fmV(λ) for λ in λs]
    V0s = [-1 / fV(λ, μs, λs .- im*c) for λ in λs]

    Id_p_U = Matrix{ComplexF64}(I, (N+1, N+1))
    for j in 1:N+1
        factor = im / (Vps[j] - Vms[j]) * fV(λs[j], μs, [λs[1:j-1]; λs[j]+1; λs[j+1:end]])
        for k in 1:N+1
            Id_p_U[j, k] += factor * (-gtheta(λs[j], λs[k], c) - gtheta(p, λs[k], c) * gtheta(s, λs[j], c))
        end
    end

    result = det(Id_p_U) / (fpV(p) - fmV(p)) / (fpV(s) - fmV(s))
    result_ln_norm = log(norm(result))
    result_angle = angle(result)

    for j in 1:N+1
        vj = (Vps[j] - Vms[j]) * V0s[j]
        result_ln_norm += log(norm(vj))
        result_angle += angle(vj)
    end

    result_ln_norm -= 0.5*log(c) + ln_gaudin_norm(ψ1) + ln_gaudin_norm(ψ2)

    return result_ln_norm, exp(im*result_angle) * im 
end

#"""
#    form_factor of the field operator <ψ1|ψ(0)|ψ2>. ψ1 and ψ2 will be normalized.
#    also correct. from T. Kojima, V. E Korepin, NA Slavnov, Comm. Math. Phys. 188, 657–689 (1997)
#"""
#function ψ0_form_factor1(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat})
#    c, L = ψ1.c, ψ1.L
#    N = length(ψ1.ns)
#
#    g(λ, μ) = im*c / (λ - μ)
#    #f(λ, μ) = im*c / (λ - μ) + 1
#    h(λ, μ) = (λ - μ + im*c) / (im*c)
#    t(λ, μ) = g(λ, μ) / h(λ, μ)
#    #a(λ) = exp(-im*L*λ / 2)
#    #d(λ) = exp( im*L*λ / 2)
#
#    μs = ψ1.quasimomenta
#    λs = ψ2.quasimomenta
#
#    detSs = zeros(ComplexF64, N+1)
#    for l in 1:N+1 
#        S = zeros(ComplexF64, (N, N))
#        js = [1:l-1 ; l+1:N+1]
#        for jx in 1:N 
#            for k in 1:N
#                j = js[jx]
#                S[jx, k] = t(μs[k], λs[j]) / h(λs[N+1], λs[j]) *
#                        prod([h(μs[ix], λs[j]) / h(λs[ix], λs[j]) for ix in 1:N]) - 
#                           t(λs[j], μs[k]) / h(λs[j], λs[N+1]) *
#                        prod([h(λs[j], μs[ix]) / h(λs[j], λs[ix]) for ix in 1:N])
#            end
#        end
#        detSs[l] = det(S)
#    end
#    Mi = sum([(-1)^(ix-1) * detSs[ix] for ix in 1:N+1])
#
#    ln_norm_FN = log(norm(Mi))
#    angle_FN = angle(Mi * exp(0.5*im*L*sum([μs; λs])) )
#
#    for m in 1:N+1
#        for j in 1:N+1
#            h_mj = h(λs[m], λs[j])
#            ln_norm_FN += log(norm(h_mj)) 
#            angle_FN += angle(h_mj)
#        end
#    end
#
#    for k in 1:N
#        for j in k+1:N+1
#            g_kj = g(λs[k], λs[j]) 
#            (j < N+1) && (g_kj *= g(μs[j], μs[k]))
#            ln_norm_FN += log(norm(g_kj))
#            angle_FN += angle(g_kj)
#        end
#    end
#
#    phase_FN = -im * exp(im*angle_FN)
#    ln_norm_FN += 0.5 * log(c)
#
#    ln_norm_FN -= ln_gaudin_norm(ψ1) + ln_gaudin_norm(ψ2) 
#    return ln_norm_FN, phase_FN
#end
