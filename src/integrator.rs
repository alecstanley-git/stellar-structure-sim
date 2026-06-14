use crate::{
    constants::{A, C, G},
    physics,
    star::Shell,
};
use std::f64::consts::PI;

// RK4 solving requires a pure function to obtain numerical derivatives without modifying the state.
// Composition is passed in as the *local* hydrogen fraction `x` for this mass shell so that the
// integrator stays correct as the core is depleted of hydrogen during evolution. Helium is taken
// as y = 1 - x - z so that the total mass fraction is always conserved.
pub fn get_derivatives(
    m: f64,
    r: f64,
    p: f64,
    l: f64,
    t: f64,
    x: f64,
    z: f64,
) -> (f64, f64, f64, f64) {
    // Always avoid absolute zero
    let safe_r = r.max(1e-10);
    let safe_m = m.max(1e-10);
    let safe_t = t.max(1e-10);
    let safe_p = p.max(1e-10);

    let y = (1.0 - x - z).max(0.0);

    // Get all the derived quantities
    let mu = physics::mean_molecular_weight(x, y, z);
    let rho = physics::calculate_density(safe_p, safe_t, mu).max(1e-10);
    let kappa = physics::calculate_opacity(rho, safe_t, x, z);
    let epsilon = physics::calculate_epsilon(rho, safe_t, x);

    // Calculate the derivatives
    let dr_dm = 1.0 / (4.0 * PI * safe_r.powi(2) * rho);
    let dp_dm = -(G * safe_m) / (4.0 * PI * safe_r.powi(4));
    let dl_dm = epsilon;

    // Radiative gradient: Nabla_rad = (3 * kappa * l * P) / (16 * pi * a * c * G * m * T^4)
    let nabla_rad = (3.0 * kappa * l * safe_p) / (16.0 * PI * A * C * G * safe_m * safe_t.powi(4));
    let nabla_ad = 0.4; // Adiabatic gradient is roughly 0.4 for monatomic ideal gas
    let nabla = nabla_rad.min(nabla_ad);

    let dt_dm = -(G * safe_m * safe_t * nabla) / (4.0 * PI * safe_r.powi(4) * safe_p);

    (dr_dm, dp_dm, dl_dm, dt_dm)
}

// Take a single RK4 step of size dm using the local hydrogen fraction `x` and metallicity `z`.
pub fn rk4_step(shell: &Shell, dm: f64, x: f64, z: f64) -> Shell {
    let (m, r, p, l, t) = (shell.m, shell.r, shell.p, shell.l, shell.t);

    // k1
    let (dr1, dp1, dl1, dt1) = get_derivatives(m, r, p, l, t, x, z);

    // k2
    let m2 = m + (0.5 * dm);
    let (dr2, dp2, dl2, dt2) = get_derivatives(
        m2,
        r + 0.5 * dm * dr1,
        p + 0.5 * dm * dp1,
        l + 0.5 * dm * dl1,
        t + 0.5 * dm * dt1,
        x,
        z,
    );

    // k3
    let (dr3, dp3, dl3, dt3) = get_derivatives(
        m2,
        r + 0.5 * dm * dr2,
        p + 0.5 * dm * dp2,
        l + 0.5 * dm * dl2,
        t + 0.5 * dm * dt2,
        x,
        z,
    );

    // k4
    let m4 = m + dm;
    let (dr4, dp4, dl4, dt4) = get_derivatives(
        m4,
        r + dm * dr3,
        p + dm * dp3,
        l + dm * dl3,
        t + dm * dt3,
        x,
        z,
    );

    let mut next_shell = Shell {
        m: m4,
        r: r + (dm / 6.0) * (dr1 + 2.0 * dr2 + 2.0 * dr3 + dr4),
        p: p + (dm / 6.0) * (dp1 + 2.0 * dp2 + 2.0 * dp3 + dp4),
        l: l + (dm / 6.0) * (dl1 + 2.0 * dl2 + 2.0 * dl3 + dl4),
        t: t + (dm / 6.0) * (dt1 + 2.0 * dt2 + 2.0 * dt3 + dt4),
        x,
        ..Default::default()
    };

    next_shell.update_derivatives(x, z);

    next_shell
}
