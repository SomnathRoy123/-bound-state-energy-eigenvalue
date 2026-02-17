include("../src/LennardJonesWKB.jl")
using .LennardJonesWKB
using PyPlot

# --- Setup Systems ---
# Gamma values derived from your previous code inputs (approximate)
hydrogen = SystemParams(21.7, "Hydrogen (H₂)")
oxygen = SystemParams(150.0, "Oxygen (O₂)")

# --- 1. Solve for Eigenvalues ---
E_H2 = find_bound_states(hydrogen)
E_O2 = find_bound_states(oxygen)

# --- 2. Visualization ---

# Plot A: Energy Levels in the Potential Well
figure(figsize=(10, 6))
x = range(0.9, stop=3.0, length=1000)
V = potential.(x)

plot(x, V, "k-", linewidth=2, label="Lennard-Jones Potential")
ylim(-1.1, 0.5)
xlim(0.8, 2.5)
xlabel("Interatomic Distance (dim-less)")
ylabel("Energy (dim-less)")
title("WKB Bound States")

# Draw Oxygen Levels
for (n, E) in enumerate(E_O2)
    # Find turning points for plotting the line width
    _, r_in, r_out = LennardJonesWKB.action_integral(E)
    plot([r_in, r_out], [E, E], "b-", alpha=0.6)
    if n % 5 == 0 # Label every 5th state to avoid clutter
        text(r_out+0.05, E, "n=$(n-1)", fontsize=8, color="blue")
    end
end

# Draw Hydrogen Levels (Red dashed)
for E in E_H2
    _, r_in, r_out = LennardJonesWKB.action_integral(E)
    plot([r_in, r_out], [E, E], "r--", linewidth=2)
end

plot([], [], "b-", label="O₂ States ($(length(E_O2)))")
plot([], [], "r--", label="H₂ States ($(length(E_H2)))")
legend()
grid(true)
savefig("energy_levels.png")

# --- 3. New Feature: WKB Wavefunctions ---
# Plotting the probability density for the n=5 state of Oxygen
if length(E_O2) > 5
    figure(figsize=(10, 6))
    E_target = E_O2[6] # n=5 (0-indexed)
    
    # Reconstruct Wavefunction
    x_wave = range(0.9, stop=1.5, length=2000)
    psi_prob = wkb_wavefunction(x_wave, E_target, oxygen)
    
    # Scale wavefunction to fit on potential plot for visualization
    scale_factor = 0.1 / maximum(psi_prob) 
    psi_shifted = (psi_prob .* scale_factor) .+ E_target
    
    plot(x, V, "k-", alpha=0.3)
    plot(x_wave, psi_shifted, "g-", linewidth=1, label="|ψ|² (WKB)")
    axhline(E_target, color="b", linestyle="--", label="Energy Level n=5")
    
    xlim(0.9, 1.4)
    ylim(E_target - 0.05, E_target + 0.15)
    title("WKB Wavefunction Reconstruction (Oxygen, n=5)")
    legend()
    savefig("wavefunction_n5.png")
end


function plot_birge_sponer(states, molecule_name)
    # Calculate energy differences (Delta E)
    delta_E = diff(states)
    n_values = 0:(length(delta_E)-1)
    
    figure()
    plot(n_values, delta_E, "o-", markersize=8, label="WKB Spacing")
    
    # Linear fit to show deviation (Anharmonicity)
    # Ideally, this curve hits 0 at the dissociation limit
    xlabel("Quantum Number n")
    ylabel("Energy Gap ΔE (E_{n+1} - E_n)")
    title("Birge-Sponer Plot: Anharmonicity of $molecule_name")
    grid(true)
    legend()
    savefig("birge_sponer_$molecule_name.png")
end

# Usage:
plot_birge_sponer(E_O2, "Oxygen")




function plot_phase_space(sys::LennardJonesWKB.SystemParams, states)
    figure(figsize=(8,8))
    x_fine = range(0.8, 2.5, length=500)
    
    for (n, E) in enumerate(states)
        # Calculate momentum p(x) for this energy
        p = [LennardJonesWKB.momentum(x, E) for x in x_fine]
        
        # Filter out real parts only (where E > V)
        mask = p .> 0
        
        # Plot top and bottom half of the orbit
        plot(x_fine[mask], p[mask], color="black", alpha=0.5)
        plot(x_fine[mask], -p[mask], color="black", alpha=0.5)
    end
    
    xlabel("Position x")
    ylabel("Momentum p")
    title("Phase Space Trajectories (Quantized Orbits)")
    grid(true)
    savefig("phase_space_orbits.png")
end

# Usage:
plot_phase_space(oxygen, E_O2)
