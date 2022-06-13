@timedtestset "test form factors" for test_ix in 1:10
    L = 11
    c = 1.21
    N = 2

    ψ0 = ground_state(c, L, N)
    ψ1 = ph_excitation(ψ0, [0], [-1])
    ψ2 = ground_state(c, L, N+1)
    ψ3 = ph_excitation(ψ2, [0], [-1])

    f = twobody_wave_function(ψ0)
    f1 = twobody_wave_function(ψ1)
    f2 = threebody_wave_function(ψ2)
    f3 = threebody_wave_function(ψ3)

    norm0 = hcubature(v -> norm(f(v...))^2, (0, 0), (L, L))[1] |> sqrt
    norm1 = hcubature(v -> norm(f1(v...))^2, (0, 0), (L, L))[1] |> sqrt
    norm2 = hcubature(v -> norm(f2(v...))^2, (0, 0, 0), (L, L, L))[1] |> sqrt
    norm3 = hcubature(v -> norm(f3(v...))^2, (0, 0, 0), (L, L, L))[1] |> sqrt

    density_form_factor_10 = 2 * quadgk(x2 -> f1(0, x2)' * f(0, x2), 0, L)[1] / (norm0 * norm1)
    ln_ρ10, phase_ρ10 = ln_ρ0_form_factor(ψ1, ψ0; p=rand())
    @test norm(exp(ln_ρ10) * phase_ρ10 - density_form_factor_10) < 1e-8 

    density_form_factor_23 = 3 * hcubature(v -> f2(0, v...)' * f3(0, v...), (0, 0), (L, L))[1] / (norm2 * norm3)
    ln_ρ23, phase_ρ23 = ln_ρ0_form_factor(ψ2, ψ3; p=rand())
    @test norm(exp(ln_ρ23) * phase_ρ23 - density_form_factor_23) < 1e-8 

    # why does it have to be sqrt(3) here?
    field_form_factor_02 = sqrt(3) * hcubature(v -> f(v...)' * f2(0, v...), (0, 0), (L, L))[1] / (norm0 * norm2) 
    ln_ψ02, phase_ψ02 = ln_ψ0_form_factor(ψ0, ψ2; s=rand(), p=rand())
    @test norm(exp(ln_ψ02) * phase_ψ02 - field_form_factor_02) < 1e-8

    field_form_factor_13 = sqrt(3) * hcubature(v -> f1(v...)' * f3(0, v...), (0, 0), (L, L))[1] / (norm1 * norm3) 
    ln_ψ13, phase_ψ13 = ln_ψ0_form_factor(ψ1, ψ3; s=rand(), p=rand())
    @test norm(exp(ln_ψ13) * phase_ψ13 - field_form_factor_13) < 1e-8

    I_form_factor_23 = 6 * quadgk(x -> f2(0, 0, x)' * f3(0, 0, x), 0, L)[1] / (norm2 * norm3)
    ln_I23, phase_I23 = ln_I0_form_factor(ψ2, ψ3; s=rand(), p=rand())
    @test norm(exp(ln_I23) * phase_I23 - I_form_factor_23) < 1e-8
end