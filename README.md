# WKB Quantization of Lennard-Jones Potential

This repository calculates the bound state energy eigenvalues of diatomic molecules (specifically $H_2$ and $O_2$) modeled by the Lennard-Jones potential using the **WKB (Wentzel–Kramers–Brillouin) Semiclassical Approximation**.

## Physics

The code solves the quantization condition:
$$\int_{x_1}^{x_2} \sqrt{2\mu(E_n - V(x))} dx = \left(n + \frac{1}{2}\right) \hbar \pi$$

Where the potential is given by the dimensionless Lennard-Jones form:
$$V(x) = 4 \left[ \left(\frac{1}{x}\right)^{12} - \left(\frac{1}{x}\right)^{6} \right]$$

## Features
* **Automatic Root Finding:** Unlike standard scripts that require manual guesses, this algorithm automatically detects the number of bound states and brackets the energy roots.
* **Adaptive Integration:** Uses `QuadGK.jl` for high-precision computation of the action integral.
* **Wavefunction Reconstruction:** Reconstructs the semiclassical probability density $|\psi(x)|^2$ in the classically allowed region.
* **Extensible:** Uses a `SystemParams` struct, making it easy to add new molecules (Nitrogen, Argon, etc.).

## Dependencies
* `QuadGK`
* `Roots`
* `PyPlot`
