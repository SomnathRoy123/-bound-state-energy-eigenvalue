module LennardJonesWKB

using QuadGK     # Standard library for adaptive integration
using Roots      # Robust root finding (add via ] add Roots)
using LinearAlgebra

export SystemParams, potential, find_bound_states, wkb_wavefunction

# --- 1. System Definition ---
struct SystemParams
    gamma::Float64  # Dimensionless quantum parameter (proportional to sqrt(mass))
    name::String
end

# The Lennard-Jones Potential: V(x) = 4(1/x^12 - 1/x^6)
potential(x) = 4 * (x^(-12) - x^(-6))

# Momentum function: p(x) = sqrt(E - V(x)) (scaled)
function momentum(x, E)
    V = potential(x)
    return E > V ? sqrt(E - V) : 0.0
end

# --- 2. Action Integral S(E) ---
# Calculates ∮ p(x) dx between turning points
function action_integral(E)
    # 1. Find turning points (roots of V(x) - E = 0)
    # For LJ, inner turning point is between 0.8 and 1.12, outer is > 1.12
    # We find roots of potential(x) - E = 0
    
    r_eq = 2^(1/6) # Equilibrium position ~1.122
    
    # Root finding wrapper
    V_diff(x) = potential(x) - E
    
    # Robustly find turning points using bracket search
    r_in = find_zero(V_diff, (0.8, r_eq))
    r_out = find_zero(V_diff, (r_eq, 10.0)) # Assuming E < 0 bound state
    
    # Integrate momentum between turning points using adaptive Gauss-Kronrod
    integral, error = quadgk(x -> sqrt(E - potential(x)), r_in, r_out)
    
    return integral, r_in, r_out
end

# --- 3. Quantization Condition ---
# Returns mismatch: S(E) - (n + 1/2)π/γ
function quantization_condition(E, n, gamma)
    integral, _, _ = action_integral(E)
    target_action = (n + 0.5) * π / gamma
    return integral - target_action
end

# --- 4. Main Solver ---
function find_bound_states(sys::SystemParams)
    states = Float64[]
    n = 0
    
    println("Solving for $(sys.name) (γ = $(sys.gamma))...")
    
    while true
        # Define the function for which we want the root: f(E) = 0
        f(E) = quantization_condition(E, n, sys.gamma)
        
        # Check if a state exists:
        # At the bottom of the well (E -> -1), the integral is 0.
        # So f(-1) = -target < 0.
        # We check if f(0) > 0 (meaning the action allows for this state before dissociation)
        if f(-0.0001) < 0
            break # No more bound states fit in the well
        end
        
        # Find root between bottom of well (-1) and dissociation limit (0)
        try
            # Search for energy E in [-0.999, -0.001]
            E_n = find_zero(f, (-0.9999, -0.0001))
            push!(states, E_n)
            println("  Found n=$n: E = $(round(E_n, digits=5))")
            n += 1
        catch e
            break
        end
    end
    return states
end

# --- 5. New Feature: Wavefunction Reconstruction ---
# Reconstructs the semiclassical probability density |ψ|²
function wkb_wavefunction(x_grid, E, sys::SystemParams)
    # Find turning points
    _, r_in, r_out = action_integral(E)
    
    psi_sq = zeros(length(x_grid))
    
    for (i, x) in enumerate(x_grid)
        if r_in < x < r_out
            p = momentum(x, E)
            # Standard WKB amplitude ∝ 1/p
            # We ignore the oscillatory phase for the probability envelope plot
            # to make it look cleaner, or we can include the cos^2 term.
            
            # Phase integral
            S_x, _ = quadgk(y -> momentum(y, E), r_in, x)
            phase = sys.gamma * S_x - π/4
            
            psi_sq[i] = (1/p) * cos(phase)^2
        else
            psi_sq[i] = 0.0 # Classical forbidden region (simple approx)
        end
    end
    
    # Normalize roughly
    psi_sq ./= sum(psi_sq) * (x_grid[2]-x_grid[1])
    return psi_sq
end

end
