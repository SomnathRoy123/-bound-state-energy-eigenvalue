include("../src/LennardJonesWKB.jl")
using .LennardJonesWKB
using PyPlot

# ============================================================
# Helper Functions (Fixing the missing definitions)
# ============================================================

# Wrapper to get turning points from the module's action_integral function
function get_turning_points(E)
    # action_integral returns (integral, r_in, r_out)
    _, r_in, r_out = LennardJonesWKB.action_integral(E)
    return (r_in, r_out)
end

# ============================================================
# Systems
# ============================================================

hydrogen = SystemParams(21.7, "Hydrogen (H₂)")
oxygen   = SystemParams(150.0, "Oxygen (O₂)")

# ============================================================
# Solve Eigenvalues
# ============================================================

println("Calculating Hydrogen states...")
E_H2 = find_bound_states(hydrogen)
tp_H2 = [get_turning_points(E) for E in E_H2]

println("Calculating Oxygen states...")
E_O2 = find_bound_states(oxygen)
tp_O2 = [get_turning_points(E) for E in E_O2]

# ============================================================
# Plot A — Energy Levels
# ============================================================

figure(figsize=(10,6))

x = range(0.85, stop=4.0, length=2000)
V = LennardJonesWKB.potential.(x)  # Added module prefix just in case

plot(x, V, "k-", linewidth=2, label="Lennard-Jones Potential")
ylim(-1.1, 0.2)
xlim(0.8, 3.5)
xlabel("Interatomic Distance (dimensionless)")
ylabel("Energy")
title("WKB Bound States")

# Oxygen (blue)
for (n, (E, (r_in, r_out))) in enumerate(zip(E_O2, tp_O2))
    plot([r_in, r_out], [E, E], "b-", alpha=0.5)
    # Label every 5th state instead of 10th for better visibility
    if n % 5 == 1 
        text(r_out + 0.05, E, "n=$(n-1)", fontsize=8, color="blue")
    end
end

# Hydrogen (red dashed)
for (E, (r_in, r_out)) in zip(E_H2, tp_H2)
    plot([r_in, r_out], [E, E], "r--", linewidth=2)
end

# Dummy plots for legend
plot(NaN, NaN, "b-", label="O₂ ($(length(E_O2)) states)")
plot(NaN, NaN, "r--", label="H₂ ($(length(E_H2)) states)")
legend()
grid(true)
savefig("plots/energy_levels.png", dpi=200) # Suggest saving to plots/ folder

# ============================================================
# WKB Wavefunction (example n=5 oxygen)
# ============================================================

if length(E_O2) > 6

    figure(figsize=(10,6))

    n = 5
    E_target = E_O2[n+1] # Julia is 1-indexed, so n=5 is index 6
    r_in, r_out = tp_O2[n+1]

    # Zoom in on the classically allowed region
    x_wave = range(r_in*0.7, r_out*1.3, length=4000)
    ψ2 = wkb_wavefunction(x_wave, E_target, oxygen)

    # Scale for visualization
    scale = 0.15 / maximum(ψ2) 
    ψplot = ψ2 .* scale .+ E_target

    plot(x, V, "k-", alpha=0.3)
    plot(x_wave, ψplot, "g-", linewidth=1.2, label="|ψ|² WKB")
    axhline(E_target, color="b", linestyle="--", label="n=$n Level")

    xlim(r_in*0.8, r_out*1.2)
    ylim(E_target - 0.1, E_target + 0.2)

    title("WKB Wavefunction Reconstruction (O₂, n=$n)")
    legend()
    savefig("plots/wavefunction_n5.png", dpi=200)
end

# ============================================================
# Birge–Sponer Plot
# ============================================================

function plot_birge_sponer(states, name)
    # Energy differences
    ΔE = diff(states)
    
    # Quantum numbers for the gaps
    v = (0:length(ΔE)-1) .+ 0.5

    figure(figsize=(7,5))
    plot(v, ΔE, "o-", markersize=6, label="WKB spacing", color="purple")

    xlabel("Quantum Number v + 1/2")
    ylabel("Energy Gap ΔE = E_{n+1} - E_n")
    title("Birge–Sponer Plot — $name")
    grid(true)
    legend()

    # Sanitized filename
    safe_name = replace(name, " " => "_")
    safe_name = replace(safe_name, "(" => "")
    safe_name = replace(safe_name, ")" => "")
    savefig("plots/birge_sponer_$(safe_name).png", dpi=200)
end

plot_birge_sponer(E_O2, "Oxygen")

# ============================================================
# Phase Space Plot
# ============================================================

function plot_phase_space(states, name)
    figure(figsize=(7,7))
    x_grid = range(0.85, 3.5, length=1500)

    # Plot a subset of states to avoid clutter
    subset_indices = 1:2:length(states) # Every 2nd state
    
    for i in subset_indices
        E = states[i]
        # Use module prefix for momentum
        p = [LennardJonesWKB.momentum(xi, E) for xi in x_grid]
        
        # Mask where momentum is real (classically allowed)
        mask = p .> 0

        plot(x_grid[mask], p[mask], "k-", alpha=0.5, linewidth=1)
        plot(x_grid[mask], -p[mask], "k-", alpha=0.5, linewidth=1)
    end

    xlabel("Position (x)")
    ylabel("Momentum (p)")
    title("Phase-Space Orbits (Quantized) — $name")
    grid(true)
    savefig("plots/phase_space_orbits.png", dpi=200)
end

plot_phase_space(E_O2, "Oxygen")
