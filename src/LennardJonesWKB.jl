module LennardJonesWKB

using QuadGK
using Roots
using LinearAlgebra

export SystemParams, potential, find_bound_states, wkb_wavefunction, harmonic_approx

# ============================================================
# 1. System Definition
# ============================================================

struct SystemParams
    gamma::Float64
    name::String
end

# Lennard-Jones potential  V(x)=4(1/x^12 − 1/x^6)
potential(x) = 4*(x^-12 - x^-6)

# momentum p(x)=sqrt(E−V)
momentum(x,E) = E > potential(x) ? sqrt(E - potential(x)) : 0.0


# ============================================================
# 2. Turning Points (ROBUST)
# ============================================================

function turning_points(E)

    f(x) = potential(x) - E
    r_eq = 2^(1/6)

    # inner root
    r_in = find_zero(f, (0.5, r_eq))

    # adaptive search for outer root
    r = r_eq*1.5
    while f(r) < 0
        r *= 1.5
        r > 1e6 && error("Outer turning point not found")
    end

    r_out = find_zero(f, (r_eq, r))
    return r_in, r_out
end


# ============================================================
# 3. Action Integral  S(E)=∫ p dx
# ============================================================

function action_integral(E)

    r_in, r_out = turning_points(E)

    integrand(x) = sqrt(max(E - potential(x), 0.0))

    S, _ = quadgk(integrand, r_in, r_out, rtol=1e-10)

    return S, r_in, r_out
end


# ============================================================
# 4. WKB Quantization
# 2γ∫pdx = π(2n+1)
# ============================================================

function quantization_condition(E, n, γ)
    S, _, _ = action_integral(E)
    return 2γ*S - π*(2n+1)
end


# ============================================================
# 5. Bound State Finder
# ============================================================

function find_bound_states(sys::SystemParams)

    println("Solving for $(sys.name)  (γ=$(sys.gamma))")

    states = Float64[]
    n = 0

    Emin = -0.9999
    Emax = -1e-8

    while true

        f(E) = quantization_condition(E, n, sys.gamma)

        # state exists only if bracketed
        if f(Emin)*f(Emax) > 0
            break
        end

        try
            E_n = find_zero(f, (Emin, Emax), Bisection())
            push!(states, E_n)
            println("n=$n   E=$(round(E_n,digits=8))")
            n += 1
        catch
            break
        end
    end

    println("Total bound states = $(length(states))")
    return states
end


# ============================================================
# 6. WKB Wavefunction (CORRECT — includes tunneling tails)
# ============================================================

function wkb_wavefunction(xgrid, E, sys::SystemParams)

    γ = sys.gamma
    S, r_in, r_out = action_integral(E)

    ψ2 = zeros(length(xgrid))

    for (i,x) in enumerate(xgrid)

        # -------- classically allowed --------
        if r_in < x < r_out

            p = sqrt(E - potential(x))

            Sx,_ = quadgk(y->sqrt(max(E - potential(y),0.0)), r_in, x)

            phase = γ*Sx - π/4
            ψ2[i] = cos(phase)^2 / p

        # -------- forbidden left --------
        elseif x < r_in

            κ = sqrt(potential(x) - E)
            Sx,_ = quadgk(y->sqrt(potential(y)-E), x, r_in)
            ψ2[i] = exp(-2γ*Sx)/κ

        # -------- forbidden right --------
        else
            κ = sqrt(potential(x) - E)
            Sx,_ = quadgk(y->sqrt(potential(y)-E), r_out, x)
            ψ2[i] = exp(-2γ*Sx)/κ
        end
    end

    # normalize
    dx = xgrid[2]-xgrid[1]
    ψ2 ./= sum(ψ2)*dx

    return ψ2
end


# ============================================================
# 7. Harmonic Approximation near minimum
# ============================================================

function harmonic_approx(n, γ)

    x0 = 2^(1/6)

    # exact curvature
    Vpp = 4*(156/x0^14 - 42/x0^8)

    ω = sqrt(Vpp)/γ

    return -1 + (n+0.5)*ω
end

end
