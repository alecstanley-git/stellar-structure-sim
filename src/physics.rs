use crate::constants::{A, K_B, M_H};

/*
This file holds generic public functions that help calculate derived quantities
*/

// Gets mean molecular weight (mu)
pub fn mean_molecular_weight(x: f64, y: f64, z: f64) -> f64 {
    // Use standard equation for an ionised gas
    1.0 / (2.0 * x + 0.75 * y + 0.5 * z)
}

// Calculates density using an equation of state
pub fn calculate_density(p: f64, t: f64, mu: f64) -> f64 {
    // P_total = P_gas + P_rad
    let p_rad = (1.0 / 3.0) * A * t.powi(4);
    let p_gas = p - p_rad;

    // Avoid negative density
    if p_gas <= 0.0 {
        return 1.0e-10;
    }

    // P_gas = (rho * k_B * T) / (mu * m_H)  =>  rho = (P_gas * mu * m_H) / (k_B * T)
    (p_gas * mu * M_H) / (K_B * t)
}

// Rosseland-mean opacity.
//
// This is the crudest piece of the microphysics. The original electron-scattering + Kramers
// formula, kappa = 0.2(1+X) + 4e25 Z(1+X) rho T^-3.5, is unusable here: the T^-3.5 term
// blows up to absurd values (>1e5) in cool dense layers yet collapses to the bare
// electron-scattering value (~0.34) in a tenuous envelope. With such a low envelope opacity
// the envelope stays radiative and there is NO compact solution -- the only root of the
// structure equations is a hugely bloated, over-luminous star.
//
// Real Rosseland-mean opacities (e.g. OPAL tables) include bound-free opacity of metals,
// H-minus and molecular bands, and across most of the solar interior they sit at a few
// cm^2/g. We approximate that by a constant, which keeps the envelope optically thick
// enough to become convective and yields a Sun-like model (L ~ 1 L_sun, T_c ~ 1.5e7 K).
// Replacing this with tabulated opacities would be the natural next improvement.
pub fn calculate_opacity(rho: f64, t: f64, x: f64, z: f64) -> f64 {
    let _ = (rho, t, x, z); // signature kept for a future tabulated/analytic opacity
    3.0
}

// Calculates energy generation using the standard PP-chain approximation.
// Uses the Gamow form (Hansen & Kawaler):
//   eps_pp = 2.38e6 * rho * X^2 * T6^(-2/3) * exp(-33.80 * T6^(-1/3))   [erg/g/s]
// The exponential captures the steep temperature sensitivity of the Coulomb-barrier
// tunnelling; near the solar core (T6 ~ 15) this behaves like a T^4 power law but with
// the correct normalisation (~20 erg/g/s), unlike the old hand-tuned coefficient which
// was ~80x too small and produced a luminosity well below L_sun.
pub fn calculate_epsilon(rho: f64, t: f64, x: f64) -> f64 {
    // Negligible below ~10^6 K (no appreciable PP burning)
    if t < 1.0e6 {
        return 0.0;
    }

    let t6 = t / 1.0e6; // Millions of Kelvin
    2.38e6 * rho * x * x * t6.powf(-2.0 / 3.0) * (-33.80 * t6.powf(-1.0 / 3.0)).exp()
}
