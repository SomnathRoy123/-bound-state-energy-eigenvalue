module LennardJonesWKB

using QuadGK
using Roots
using LinearAlgebra

export SystemParams, potential, momentum, find_bound_states, wkb_wavefunction, action_integral

# ============================================================
# 1. System Definition
# ============================================================

struct SystemParams
    gamma::Float64
    name::String
end

# Lennard-Jones potential  V(x)=4(1/x^12 − 1/x^6)
potential(x) = 4*(x^-12 - x^-6)

# Momentum p(x)=sqrt(E−V)
function momentum(x, E)
    V = potential(x)
    return E > V ? sqrt(E - V) : 0.0
end

# ============================================================
# 2. Turning Points (ROBUST)
# ============================================================

function turning_points(E)
    f(x) = potential(x) - E
    r_eq = 2^(1/6)

    # Inner root (must be between 0.8 and r_eq)
    r_in = find_zero(f, (0.8, r_eq))

    # Adaptive search for outer root
    r_search = r_eq * 1.1
    while f(r_search) < 0
        r_search *= 1.2
        if r_search > 50.0 
            error("Outer turning point not found for E=$E")
        end
    end

    r_out = find_zero(f, (r_eq, r_search))
    return r_in, r_out
end

# ============================================================
# 3. Action Integral  S(E) = ∫ p dx
# ============================================================

function action_integral(E)
    r_in, r_out = turning_points(E)
    
    # Integrand is momentum p(x)
    integrand(x) = sqrt(max(E - potential(x), 0.0))
    
    S, _ = quadgk(integrand, r_in, r_out, rtol=1e-8)
    return S, r_in, r_out
end

# ============================================================
# 4. WKB Quantization Condition
# 2γ ∫ p dx = π(2n+1) -> Return mismatch
# ============================================================

function quantization_mismatch(E, n, gamma)
    S, _, _ = action_integral(E)
    # Bohr-Sommerfeld Condition: 2*gamma*S = (n + 0.5)*2π? 
    # Actually standard WKB is: integral p dx = (n + 0.5) * pi / gamma
    # Let's use the form: 2*gamma*S - (2n+1)*pi = 0
    return 2 * gamma * S - (2*n + 1) * π
end

# ============================================================
# 5. Bound State Finder
# ============================================================

function find_bound_states(sys::SystemParams)
    println("Solving for $(sys.name) (γ=$(sys.gamma))...")
    
    states = Float64[]
    n = 0
    
    # We search from bottom of well (-1.0) up to dissociation (0.0)
    E_min = -0.9999
    E_max = -0.0001
    
    while true
        # Define function f(E) that is 0 when E is the n-th eigenvalue
        f(E) = quantization_mismatch(E, n, sys.gamma)
        
        # Check if the n-th state even exists in the potential well
        # At E_min (bottom), action is ~0, so mismatch is negative.
        # At E_max (top), action is large.
        if f(E_min) * f(E_max) > 0
            # If signs are same, likely no root for this n (or n is too high)
            break
        end

        try
            E_n = find_zero(f, (E_min, E_max), Bisection())
            push!(states, E_n)
            println("  Found n=$n: E=$(round(E_n, digits=6))")
            n += 1
            # Optimization: Next state must be higher energy
            E_min = E_n + 1e-6 
        catch
            break
        end
    end

    println("Total bound states found: $(length(states))")
    return states
end

# ============================================================
# 6. WKB Wavefunction (Reconstruction)
# ============================================================

function wkb_wavefunction(x_grid, E, sys::SystemParams)
    gamma = sys.gamma
    _, r_in, r_out = action_integral(E)
    
    psi_sq = zeros(length(x_grid))
    epsilon = 1e-3 # Softening parameter to avoid division by zero at turning points

    for (i, x) in enumerate(x_grid)
        if x > r_in && x < r_out
            # Classically Allowed Region
            p = sqrt(E - potential(x))
            
            # Compute phase integral from r_in to x
            S_x, _ = quadgk(y -> sqrt(E - potential(y)), r_in, x, rtol=1e-5)
            
            phase = gamma * S_x - π/4
            psi_sq[i] = (cos(phase)^2) / (p + epsilon)
            
        elseif x <= r_in
            # Forbidden Left (Exponential Decay)
            kappa = sqrt(potential(x) - E)
            dS, _ = quadgk(y -> sqrt(potential(y) - E), x, r_in, rtol=1e-5)
            psi_sq[i] = (exp(-2 * gamma * dS)) / (kappa + epsilon)
            
        else
            # Forbidden Right (Exponential Decay)
            kappa = sqrt(potential(x) - E)
            dS, _ = quadgk(y -> sqrt(potential(y) - E), r_out, x, rtol=1e-5)
            psi_sq[i] = (exp(-2 * gamma * dS)) / (kappa + epsilon)
        end
    end

    # Normalize
    dx = x_grid[2] - x_grid[1]
    norm = sum(psi_sq) * dx
    return psi_sq ./ norm
end

end
