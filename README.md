# Semiclassical Bound States of the Lennard-Jones Potential

![Language](https://img.shields.io/badge/Language-Julia_1.9+-purple.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Status](https://img.shields.io/badge/Status-Research_Prototype-orange.svg)

This repository computes vibrational bound states of a diatomic molecule using the **WKB (Wentzel–Kramers–Brillouin) semiclassical approximation**.

Instead of solving the Schrödinger equation on a spatial grid, the code determines quantum energy levels from **classical phase-space quantization**. The implementation is written in **Julia** and serves as a minimal computational molecular spectroscopy model.

## ⚛️ Physical Model

We study a particle moving in the reduced Lennard-Jones potential:

$$
V(x) = 4\left(x^{-12} - x^{-6}\right)
$$

**Key Properties:**
- Strong repulsive core at short distance.
- Attractive van-der-Waals tail $V(x)\sim -1/x^6$.
- Finite number of bound vibrational states.

The dimensionless Schrödinger equation becomes:

$$
-\frac{1}{\gamma^2}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
$$

The parameter $\gamma$ controls the classicality of the molecule:

$$
\gamma = \frac{\sqrt{2\mu D_e}\,r_0}{\hbar}
$$

| Regime | Interpretation |
| :--- | :--- |
| $\gamma \ll 1$ | Strongly quantum |
| $\gamma \sim 10$ | Moderately quantum (e.g., $H_2$, $O_2$) |
| $\gamma \gg 1$ | Semiclassical molecular limit |

---

## 📐 WKB Quantization Method

We define the classical momentum as:

$$
p(x) = \sqrt{E - V(x)}
$$

The turning points $x_1, x_2$ satisfy $V(x) = E$. The **Bohr–Sommerfeld quantization rule** gives the discrete vibrational energies:

$$
2\gamma \int_{x_1}^{x_2} p(x)\,dx = \pi(2n+1)
$$

The program numerically solves this integral equation to obtain all bound states $E_n$.

---

## 📊 Physics Demonstrated

### 1. Harmonic Region (Low Energy)
Near the potential minimum, the spectrum is equally spaced:

$$
E_n \approx -1 + \frac{\omega}{\gamma}\left(n+\frac12\right)
$$

### 2. Anharmonic Vibrational Spectrum
As energy increases, the potential widens, causing the energy spacing to decrease (Birge–Sponer behavior):

$$
\Delta E_n = E_{n+1} - E_n
$$

### 3. Near Dissociation (LeRoy–Bernstein Law)
For a long-range $1/x^6$ tail, the levels cluster near the dissociation limit ($E=0$) following:

$$
E_n \propto -(n_D-n)^3
$$

### 4. Semiclassical Interpretation
Each bound state corresponds to a quantized phase-space orbit where the area is an integer multiple of Planck's constant:

$$
\oint p\,dx = (n+\tfrac12)h
$$



## 📂 Repository Structure

```text
WKB-LennardJones/
├── src/
│   └── LennardJonesWKB.jl    # Core WKB solver implementation
├── scripts/
│   └── run_analysis.jl       # Analysis and plotting script
├── test/                     # Unit tests for correctness
└── README.md                 # Project documentation

```

## 🚀 Installation & Usage

1. **Install Julia** (≥ 1.9)
2. **Clone the repository** and instantiate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
## Run the Analysis:

```julia
julia scripts/run_analysis.jl
```
## 👤 Author

**Somnath Roy**
* **Institution:** IISER Kolkata (BS-MS Physics)
* **Research Interests:** Dynamical Systems, Active Matter, Statistical Physics
* **GitHub:** [@SomnathRoy123](https://github.com/SomnathRoy123)

---

> **Scientific Purpose:** This project demonstrates how molecular vibrational spectra emerge from classical mechanics through semiclassical quantization. It serves as a pedagogical tool for connecting classical trajectories to quantum states.
