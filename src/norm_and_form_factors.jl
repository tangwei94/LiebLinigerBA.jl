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
    target = :density or :current
"""
function ln_ρ0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0, target=:density)

    if ψ1.ns == ψ2.ns
        if target == :density
            return log(length(ψ1.ns) / ψ1.L), 1.0
        elseif target == :current 
            return log(abs(momentum(ψ1)) / ψ1.L), sign(momentum(ψ1))
        end 
    end

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
    V0s = [fV(λ, λs .+ im*c, μs) for λ in λs]

    Id_p_U = Matrix{ComplexF64}(I, (N+1, N+1))
    for j in 1:N
        factor = im * (μs[j] - λs[j]) / (Vps[j] - Vms[j]) * fV(λs[j], μs, [λs[1:j-1]; μs[j]; λs[j+1:end]]) 
        for k in 1:N
            Id_p_U[j, k] += factor * (-gtheta(λs[j], λs[k], c) + gtheta(p, λs[k], c))
        end
    end

    result = det(Id_p_U) / (fpV(p) - fmV(p)) 
    result_ln_norm = log(norm(result))
    result_angle = angle(result)

    for j in 1:N
        vj = (Vps[j] - Vms[j]) * V0s[j]
        result_ln_norm += log(norm(vj))
        result_angle += angle(vj)
    end

    if target == :density
        result_ln_norm += log(norm(sum(μs) - sum(λs)))
        additional_sign = sign(sum(μs) - sum(λs))
    elseif target == :current 
        result_ln_norm += log(norm(sum(μs .^ 2) - sum(λs .^ 2)))
        additional_sign = sign(sum(μs .^ 2) - sum(λs .^ 2))
    end

    result_ln_norm -= ln_gaudin_norm(ψ1) + ln_gaudin_norm(ψ2)

    # here I have to add an additional phase factor -im, otherwise the result seems incorrect
    return result_ln_norm, exp(im*result_angle) * additional_sign * (-im) 
end

"""
    kacmoody(q::Integer, sector::Symbol, v::Real)

    microscopic realization of the kacmoody generator j_{q}.
    sector = :holomorphic, :antiholomorphic, :none
"""
function kacmoody(ϕ1::LLBAState{<:AbstractFloat}, ϕ2::LLBAState{<:AbstractFloat}, sector::Symbol, v::Real)
    if length(ϕ1.ns) != length(ϕ2.ns)
        return 0
    end
    if ϕ1.ns == ϕ2.ns
        return 0
    end

    sgn = Float64((sector==:holomorphic) - (sector==:antiholomorphic))
    lnnorm_ρ, phase_ρ = ln_ρ0_form_factor(ϕ1, ϕ2; target=:density)
    lnnorm_j, phase_j = ln_ρ0_form_factor(ϕ1, ϕ2; target=:current)
    ρ = exp(lnnorm_ρ) * phase_ρ 
    j = exp(lnnorm_j) * phase_j 
    return ϕ1.L * (v * ρ + sgn * j)  
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
function ln_ψ0_form_factor(ψ1::LLBAState{<:AbstractFloat}, ψ2::LLBAState{<:AbstractFloat}; p::Real=0, s::Real=0)
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

    # here I have to add an additional phase factor -im, otherwise the result seems incorrect
    return result_ln_norm, exp(im*result_angle) * im * (-im) 
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