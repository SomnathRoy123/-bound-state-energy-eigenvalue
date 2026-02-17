# Theory: Semiclassical Bound States of the Lennard-Jones Potential

This project computes vibrational bound states of a diatomic molecule modeled by the Lennard-Jones (12-6) potential using the semiclassical WKB approximation.

---

## 1. Lennard-Jones Potential

We work in reduced units

\[
V(x) = 4\left(x^{-12} - x^{-6}\right)
\]

Properties:

- Minimum at \( x_0 = 2^{1/6} \)
- Well depth \( V_{\min} = -1 \)
- Long-range tail \( V(x) \sim -4/x^6 \)

The \(1/x^6\) attraction determines the dissociation-limit spectrum.

---

## 2. Dimensionless Schrödinger Equation

The radial vibrational equation can be written in scaled form

\[
-\frac{1}{\gamma^2}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
\]

where

\[
\gamma = \frac{\sqrt{2\mu D_e}\, r_0}{\hbar}
\]

is the semiclassical parameter.

Interpretation:

| Regime | Meaning |
|------|------|
| \( \gamma \ll 1 \) | Strongly quantum |
| \( \gamma \gg 1 \) | Semiclassical / molecular limit |

Heavy molecules → large γ → WKB valid.

---

## 3. WKB Quantization Condition

Define classical momentum

\[
p(x)=\sqrt{E-V(x)}
\]

Turning points satisfy

\[
V(x)=E
\]

The Bohr-Sommerfeld rule gives

\[
2\gamma \int_{x_1}^{x_2} p(x)\,dx = \pi(2n+1)
\]

This determines discrete vibrational energies \(E_n\).

---

## 4. Harmonic Limit Near the Minimum

Expanding the potential around equilibrium:

\[
V(x) \approx -1 + \frac{1}{2}V''(x_0)(x-x_0)^2
\]

The spectrum becomes

\[
E_n \approx -1 + \frac{\omega}{\gamma}\left(n+\frac{1}{2}\right)
\]

valid only for low quantum numbers.

---

## 5. Near-Dissociation Spectrum (LeRoy–Bernstein Law)

For potentials

\[
V(x)\sim -\frac{C_6}{x^6}
\]

the high-lying vibrational levels obey

\[
E_n \propto -(n_D-n)^3
\]

Consequences:

* Level spacing decreases rapidly near dissociation
* Birge–Sponer plot becomes nonlinear
* Spectrum strongly anharmonic

This behavior is a diagnostic for correct semiclassical implementation.

---

## 6. Physical Interpretation

The WKB approximation treats vibrational states as quantized classical oscillations:

\[
\oint p\,dx = \text{integer} \times h
\]

Thus bound states correspond to quantized phase-space area rather than wavefunction nodes alone.

The number of bound states scales as

\[
N \propto \gamma
\]

which matches molecular spectroscopy expectations.
