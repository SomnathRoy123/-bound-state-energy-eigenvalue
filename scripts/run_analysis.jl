include("../src/LennardJonesWKB.jl")
using .LennardJonesWKB
using PyPlot

# ============================================================
# Systems
# ============================================================

hydrogen = SystemParams(21.7, "Hydrogen (H₂)")
oxygen   = SystemParams(150.0, "Oxygen (O₂)")

# ============================================================
# Solve Eigenvalues
# ============================================================

E_H2 = find_bound_states(hydrogen)
E_O2 = find_bound_states(oxygen)

tp_H2 = [turning_points(E) for E in E_H2]
tp_O2 = [turning_points(E) for E in E_O2]

# ============================================================
# Plot A — Energy Levels
# ============================================================

figure(figsize=(10,6))

x = range(0.85, stop=4.0, length=2000)
V = potential.(x)

plot(x, V, "k-", linewidth=2, label="Lennard-Jones Potential")
ylim(-1.1, 0.2)
xlim(0.8, 3.5)
xlabel("Interatomic Distance (dimensionless)")
ylabel("Energy")
title("WKB Bound States")

# Oxygen (blue)
for (n,(E,(r_in,r_out))) in enumerate(zip(E_O2,tp_O2))
    plot([r_in,r_out],[E,E],"b-",alpha=0.5)
    if n % 10 == 1
        text(r_out+0.05,E,"n=$(n-1)",fontsize=8,color="blue")
    end
end

# Hydrogen (red dashed)
for (E,(r_in,r_out)) in zip(E_H2,tp_H2)
    plot([r_in,r_out],[E,E],"r--",linewidth=2)
end

plot(NaN,NaN,"b-",label="O₂ ($(length(E_O2)) states)")
plot(NaN,NaN,"r--",label="H₂ ($(length(E_H2)) states)")
legend()
grid(true)
savefig("energy_levels.png",dpi=200)

# ============================================================
# WKB Wavefunction (example n=5 oxygen)
# ============================================================

if length(E_O2) > 6

    figure(figsize=(10,6))

    n = 5
    E_target = E_O2[n+1]
    r_in,r_out = tp_O2[n+1]

    x_wave = range(r_in*0.7, r_out*1.3, length=4000)
    ψ2 = wkb_wavefunction(x_wave,E_target,oxygen)

    scale = 0.1/maximum(ψ2)
    ψplot = ψ2 .* scale .+ E_target

    plot(x,V,"k-",alpha=0.3)
    plot(x_wave,ψplot,"g-",linewidth=1.2,label="|ψ|² WKB")
    axhline(E_target,color="b",linestyle="--",label="n=5")

    xlim(r_in*0.9,r_out*1.1)
    ylim(E_target-0.08,E_target+0.15)

    title("WKB Wavefunction (O₂, n=5)")
    legend()
    savefig("wavefunction_n5.png",dpi=200)
end

# ============================================================
# Birge–Sponer Plot
# ============================================================

function plot_birge_sponer(states,name)

    ΔE = diff(states)
    v = (0:length(ΔE)-1) .+ 0.5

    figure(figsize=(7,5))
    plot(v,ΔE,"o-",markersize=6,label="WKB spacing")

    xlabel("v = n + 1/2")
    ylabel("ΔEₙ = Eₙ₊₁ − Eₙ")
    title("Birge–Sponer Plot — $name")
    grid(true)
    legend()

    savefig("birge_sponer_$(replace(name,' '=>'_')).png",dpi=200)
end

plot_birge_sponer(E_O2,"Oxygen")

# ============================================================
# Phase Space Plot
# ============================================================

function plot_phase_space(states,name)

    figure(figsize=(7,7))
    x = range(0.85,3.5,length=1500)

    for (n,E) in enumerate(states[1:10:end]) # avoid clutter
        p = [momentum(xi,E) for xi in x]
        mask = p .> 0

        plot(x[mask], p[mask], color="black", alpha=0.4)
        plot(x[mask],-p[mask], color="black", alpha=0.4)
    end

    xlabel("x")
    ylabel("p(x)")
    title("Semiclassical Phase-Space Orbits — $name")
    grid(true)
    savefig("phase_space_orbits.png",dpi=200)
end

plot_phase_space(E_O2,"Oxygen")
