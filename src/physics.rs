use crate::constants::{A, C, G, K_B, M_H};
use std::f64::consts::PI;

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
use crate::opacity::OpacityTable;

pub fn calculate_opacity(rho: f64, t: f64, x: f64, z: f64, table: &OpacityTable, aesopus_blend: f64) -> f64 {
    // Electron scattering
    let kappa_es = 0.2 * (1.0 + x);
    
    // Free-free opacity
    let kappa_ff = 1.75e22 * (1.0 + x) * (1.0 - z) * rho * t.powf(-3.5);
    
    // Bound-free opacity (Kramers) with guillotine
    let t6 = t / 1.0e6;
    let guillotine = 1.0 + (t6 * t6);
    let kappa_bf = 4.3e25 * z * (1.0 + x) * rho * t.powf(-3.5) / guillotine;
    
    let kappa_k = kappa_ff + kappa_bf;
    let kappa_es = 0.2 * (1.0 + x);
    let kappa_interior_raw = kappa_k + kappa_es;

    // H-minus opacity (used ONLY as a limiter for the Kramers formula to prevent low-T explosions)
    let kappa_h_minus = 2.5e-32 * (z / 0.02) * rho.sqrt() * t.powi(9);
    let old_kappa = (kappa_h_minus * kappa_interior_raw) / (kappa_h_minus + kappa_interior_raw).max(1e-10);

    let log_t = t.log10();
    
    let kappa_low_t = if log_t < 4.5 {
        let t_6_table = t / 1e6;
        let log_r = (rho / t_6_table.powi(3)).log10();
        let log_kappa = table.get_log_kappa(x, log_t, log_r);
        10.0f64.powf(log_kappa)
    } else {
        0.0 
    };

    let kappa = if log_t <= 4.1 {
        kappa_low_t
    } else if log_t >= 4.4 {
        old_kappa
    } else {
        // Smoothstep blending from 4.1 to 4.4
        let w = (log_t - 4.1) / 0.3;
        let blend = w * w * (3.0 - 2.0 * w);
        let log_k_low = kappa_low_t.log10();
        let log_k_high = old_kappa.log10();
        10.0f64.powf(log_k_low * (1.0 - blend) + log_k_high * blend)
    };

    // Global relaxation blend to ease the solver into the new physics
    let final_kappa = if aesopus_blend < 1e-5 {
        old_kappa
    } else if aesopus_blend > 0.99999 {
        kappa
    } else {
        let log_old = old_kappa.log10();
        let log_new = kappa.log10();
        10.0f64.powf(log_old * (1.0 - aesopus_blend) + log_new * aesopus_blend)
    };

    final_kappa.max(1e-4)
}

// Calculates energy generation using the standard PP-chain approximation.
// Uses the Gamow form (Hansen & Kawaler):
//   eps_pp = 2.38e6 * rho * X^2 * T6^(-2/3) * exp(-33.80 * T6^(-1/3))   [erg/g/s]
// The exponential captures the steep temperature sensitivity of the Coulomb-barrier
// tunnelling; near the solar core (T6 ~ 15) this behaves like a T^4 power law but with
// the correct normalisation (~20 erg/g/s), unlike the old hand-tuned coefficient which
// was ~80x too small and produced a luminosity well below L_sun.
pub fn calculate_epsilon(rho: f64, t: f64, x: f64, z: f64) -> f64 {
    // Negligible below ~10^6 K (no appreciable PP burning)
    if t < 1.0e6 {
        return 0.0;
    }

    let t6 = t / 1.0e6; // Millions of Kelvin
    
    // PP Chain
    let eps_pp = 3.3e6 * rho * x * x * t6.powf(-2.0 / 3.0) * (-33.80 * t6.powf(-1.0 / 3.0)).exp();
    
    // CNO Cycle
    let x_cno = 0.7 * z;
    let eps_cno = 8.67e27 * rho * x * x_cno * t6.powf(-2.0 / 3.0) * (-152.28 * t6.powf(-1.0 / 3.0)).exp();
    
    eps_pp + eps_cno
}

// Calculates gravitational energy generation rate epsilon_grav.
// Using the thermodynamic relation: eps_grav = -T ds/dt = c_P T (nabla_ad/P dP/dt - 1/T dT/dt)
// assuming an ideal monatomic gas (nabla_ad = 0.4, c_P = 5/2 k_B / mu m_H).
pub fn calculate_epsilon_grav(t: f64, t_old: f64, p: f64, p_old: f64, mu: f64, dt_sec: f64) -> f64 {
    if dt_sec <= 0.0 {
        return 0.0;
    }
    let c_p = 2.5 * K_B / (mu * M_H);
    let nabla_ad = 0.4;
    let dp_dt = (p - p_old) / dt_sec;
    let dt_dt = (t - t_old) / dt_sec;
    c_p * t * ((nabla_ad / p) * dp_dt - (1.0 / t) * dt_dt)
}

// Calculates the true temperature gradient (nabla) using Mixing Length Theory (MLT)
pub fn calculate_nabla(t: f64, p: f64, m: f64, r: f64, l: f64, rho: f64, kappa: f64, mu: f64, alpha_mlt: f64) -> f64 {
    let nabla_rad = (3.0 * kappa * l * p) / (16.0 * PI * A * C * G * m * t.powi(4));
    let nabla_ad = 0.4;
    
    if nabla_rad <= nabla_ad {
        return nabla_rad;
    }
    
    // Convection is active: solve MLT cubic
    // g = G * m / r^2
    let g = (G * m) / r.powi(2).max(1e-10);
    
    // Pressure scale height H_P = P / (rho * g)
    let hp = p / (rho * g).max(1e-10);
    
    // Mixing length l_m = alpha * H_P
    let lm = alpha_mlt * hp;
    
    // Specific heat c_P for ideal monatomic gas
    let cp = 2.5 * K_B / (mu * M_H);
    
    // Dimensionless U parameter
    let u = (3.0 * A * C * t.powi(3)) / (cp * rho.powi(2) * kappa * lm.powi(2)) * (8.0 * hp / g.max(1e-10)).sqrt();
    
    let w = nabla_rad - nabla_ad;
    
    // We solve the standard cubic for the convective efficiency root x:
    // (9/8) U x^3 + x^2 + 2Ux - W = 0
    // We use Newton-Raphson starting from a reasonable guess.
    let mut x = w.sqrt().min(1.0); // Safe initial guess
    for _ in 0..20 {
        let f = (9.0 / 8.0) * u * x.powi(3) + x.powi(2) + 2.0 * u * x - w;
        let df = (27.0 / 8.0) * u * x.powi(2) + 2.0 * x + 2.0 * u;
        let dx = f / df.max(1e-30);
        x -= dx;
        if dx.abs() < 1e-6 {
            break;
        }
    }
    
    // Constrain x to positive values to avoid numerical explosions
    x = x.max(0.0);
    
    // True nabla
    nabla_ad + x.powi(2) + 2.0 * u * x
}
