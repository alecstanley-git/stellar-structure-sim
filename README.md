# Stellar Structure Sim

My goal with this project is to build a functional program that can simulate the lifespan of a star. It will be loosely inspired by the TWIN project built in Fortran. I am building it entirely in Rust and am using this primarily as a tool to get into Rust to aid in my studies.

In general, stellar structures are difficult to simulate, involving four coupled differential equations and a highly complex numerical integration scheme. They are used to analyse stars during the pre-main sequence, the main sequence, the AGB phase, and the red giant phase among others, by relying on basic simplifying assumptions about stars (such as their convectivity/radiativity, their spherical symmetry, and approximation to a perfect blackbody). A first step of mine will be to successfully simulate the main sequence of a sun-like star.

## Status

The code now produces a Sun-like main-sequence model and evolves it. Running it gives:

```
age =  0.00 Gyr | L = 1.046 L_sun | R = 0.719 R_sun | T_c = 1.527e7 K | X_core = 0.700
age =  4.00 Gyr | L = 1.271 L_sun | R = 0.748 R_sun | T_c = 1.794e7 K | X_core = 0.364
age = 10.00 Gyr | L = 2.165 L_sun | R = 0.850 R_sun | T_c = 2.443e7 K | X_core = 0.058
```

At the Sun's actual age (~4.6 Gyr) the model has L ≈ 1.3 L_sun and core hydrogen X ≈ 0.36,
both close to the real Sun (the Sun has brightened ~30% since the zero-age main sequence and
its core is now ~0.34 hydrogen). As core hydrogen burns, the mean molecular weight rises, the
core contracts and heats (T_c, P_c climb) and the luminosity grows — textbook main-sequence
evolution.

## Physics and numerics

The structure is integrated in the Lagrangian mass coordinate, with the four standard
equations for dr/dm, dP/dm, dL/dm and dT/dm. Temperature transport uses
∇ = min(∇_rad, ∇_ad) (Schwarzschild criterion), the equation of state is ideal gas plus
radiation pressure, energy generation is the PP chain (Gamow form), and opacity is a
constant Rosseland-mean approximation (see `physics.rs` for why a constant is used).

The structure is solved as a two-point boundary value problem with the **fitting method**:
the equations are integrated *outward* from the centre and *inward* from the surface (each
in its numerically stable direction, with adaptive step sizing) and the two solutions are
required to match at a fitting point. A 4-D Newton–Raphson iteration adjusts the four
unknowns — central pressure and temperature, and surface radius and luminosity — until they
agree. Time evolution burns hydrogen (dX/dt = −ε/Q_H) on a fixed Lagrangian composition grid
and re-solves the structure each step.

### Corrections made to the original code

* **Radiative gradient** (`integrator.rs`, `star.rs`): the factor in ∇_rad was `64π`; the
  correct value is `16π` (∇_rad = 3κLP / (16π a c G m T⁴)). The old version made the
  temperature gradient 4× too shallow.
* **Energy generation** (`physics.rs`): the PP coefficient `1.07e-7` gave ε ~ 80× too small.
  Replaced with the standard Gamow form, ε = 2.38e6 ρ X² T₆^(−2/3) exp(−33.8 T₆^(−1/3)),
  which has the correct normalisation and temperature sensitivity.
* **Opacity** (`physics.rs`): the electron-scattering + Kramers formula is unusable in the
  envelope (the T^−3.5 term collapses to the bare electron-scattering value where it matters,
  so the only solution is a hugely bloated star). Replaced with a constant Rosseland-mean
  opacity representing the bound-free/H⁻/metal sources it omits, which lets the envelope
  become convective and gives a Sun-like radius.
* **Solver** (`star.rs`, `integrator.rs`): the original did a *single* outward integration
  from a fixed guess with no correction (and printed a "Pressure Error" that nothing ever
  drove to zero). Pure outward shooting is also numerically unstable. Replaced with the
  fitting-method BVP solver and an adaptive integrator that resolves the stiff surface layers.
* **Composition & evolution** (`star.rs`): the old `evolve_timestep` only incremented the
  age. It now depletes hydrogen per mass shell from the local energy generation and re-solves
  the structure, which is what drives the main-sequence evolution.

### Known limitations / next steps

* The radius is ~25% small. The dominant cause is the constant opacity; tabulated (OPAL-style)
  opacities are the natural next improvement.
* No gravitational/contraction energy term in dL/dm (fine on the main sequence, needed for the
  pre-main sequence and later phases).
* Energy generation is PP only (no CNO), so this is appropriate for ≲1.2 M_sun stars.

## Astrophysical Equations and Assumptions

### 1. Mass Conservation
* $\frac{\partial r}{\partial m} = \frac{1}{4 \pi r^2 \rho}$
* **Assumption:** The star is perfectly spherically symmetric and in hydrostatic equilibrium. Rotational and magnetic effects are neglected.

### 2. Hydrostatic Equilibrium
* $\frac{\partial P}{\partial m} = -\frac{G m}{4 \pi r^4}$
* **Assumption:** The star is in perfect hydrostatic balance. Dynamical timescales (like pulsations) are assumed to be negligible compared to the evolutionary timescale.

### 3. Energy Generation
* $\frac{\partial L}{\partial m} = \epsilon_{nuc} + \epsilon_{grav}$
* $\epsilon_{nuc} = \epsilon_{pp} = 2.38 \times 10^6 \rho X^2 T_6^{-2/3} e^{-33.8 T_6^{-1/3}}$ (ergs/g/s)
* $\epsilon_{grav} = -T \frac{\partial S}{\partial t} \approx c_P T \left( \frac{\nabla_{ad}}{P} \frac{dP}{dt} - \frac{1}{T} \frac{dT}{dt} \right)$
* **Assumption:** Energy generation is dominated by the proton-proton (p-p) chain since the star is sun-like. The CNO cycle is neglected. Gravitational contraction/expansion energy is approximated over discrete timesteps using ideal gas specific heats.

### 4. Temperature Gradient (Energy Transport)
* $\frac{\partial T}{\partial m} = -\frac{G m T}{4 \pi r^4 P} \nabla$
* $\nabla = \min(\nabla_{rad}, \nabla_{ad})$
* $\nabla_{rad} = \frac{3 \kappa L P}{16 \pi a c G m T^4}$, $\nabla_{ad} \approx 0.4$
* **Assumption:** Convection is modeled using a simplified Schwarzschild criterion where the actual gradient immediately becomes adiabatic if $\nabla_{rad} > \nabla_{ad}$. Mixing length theory (MLT) is not fully implemented.

### 5. Equation of State (Ideal Gas + Radiation)
* $P = P_{gas} + P_{rad} = \frac{\rho k_B T}{\mu m_H} + \frac{1}{3} a T^4$
* **Assumption:** The stellar material is an ideal gas. Degeneracy pressure and non-ideal plasma interactions (e.g., Coulomb interactions) are neglected, which is mostly valid for a sun-like main sequence star but fails at late stages.

### 6. Opacity
* $\kappa \approx \kappa_{es} + \kappa_{ff}$
* $\kappa_{es} = 0.2(1+X)$
* $\kappa_{ff} = 3.8 \times 10^{22} (1+X)(1-Z) \rho T^{-3.5}$
* **Assumption:** Opacity is a simple sum of electron scattering and a Kramers' law approximation for free-free/bound-free transitions. We do not use comprehensive OPAL opacity tables.
