# Numerical Methods

This project solves the semiclassical quantization condition numerically rather than integrating the Schrödinger equation directly.

---

## 1. Turning Point Bracketing

Turning points satisfy

\[
V(x) - E = 0
\]

The Lennard-Jones potential creates two roots:

* inner turning point near the repulsive wall
* outer turning point diverging as \(E \to 0^-\)

Because the outer root moves far outward near dissociation, a fixed bracket fails.

We therefore use **adaptive bracketing**:

1. Start near equilibrium distance
2. Expand search interval exponentially
3. Locate root via bisection

This guarantees convergence for all bound states.

---

## 2. Action Integral Evaluation

The action integral is

\[
S(E) = \int_{x_1}^{x_2}\sqrt{E-V(x)}\,dx
\]

The integrand has square-root singularities:

\[
p(x)\sim \sqrt{x-x_t}
\]

Standard quadrature performs poorly here.

We use **Gauss–Kronrod adaptive quadrature (QuadGK)** because it:

* handles endpoint singularities
* provides error estimation
* adapts sampling density automatically

---

## 3. Root Finding for Energy Levels

Bound states satisfy

\[
F(E) = 2\gamma S(E) - \pi(2n+1)=0
\]

We solve using bracketing methods (Bisection):

Advantages:
- Guaranteed convergence
- No derivative required
- Robust near dissociation limit

---

## 4. WKB Wavefunction Reconstruction

Classically allowed region:

\[
|\psi|^2 \propto \frac{1}{p(x)}\cos^2\left(\gamma\int p\,dx - \pi/4\right)
\]

Forbidden region:

\[
|\psi|^2 \propto \frac{1}{\kappa(x)}\exp\!\left[-2\gamma\int \kappa(x)\,dx\right]
\]

where

\[
\kappa(x)=\sqrt{V(x)-E}
\]

This produces physically correct tunneling tails.

---

## 5. Computational Complexity

Let \(N\) be number of bound states.

Each level requires:
- root finding
- action integral evaluation

Thus cost scales approximately

\[
O(N_{\text{states}})
\]

which is far cheaper than grid Schrödinger solvers \(O(N_{grid}^2)\).

This makes WKB ideal for molecular spectroscopy regimes with many levels.
