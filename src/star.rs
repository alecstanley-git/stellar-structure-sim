use std::f64::consts::PI;

use crate::{
    constants::{G, M_SUN, Q_H, SEC_PER_YEAR, SIGMA},
    integrator::{get_derivatives, rk4_step},
    parameters::{COMP_GRID, STEP_TOL},
    physics::{calculate_density, calculate_epsilon, calculate_opacity, mean_molecular_weight},
};

#[derive(Debug, Clone, Copy, Default)]
pub struct Shell {
    // Independent variable
    pub m: f64, // Enclosed mass

    // Dependent variables
    pub r: f64, // Radius
    pub p: f64, // Pressure
    pub l: f64, // Luminosity
    pub t: f64, // Temp

    // Local composition (hydrogen mass fraction at this shell)
    pub x: f64,

    // Derived quantities
    pub rho: f64,     // Density
    pub kappa: f64,   // Opacity
    pub epsilon: f64, // Energy generation rate
    pub nabla: f64,   // Temperature gradient

    // Derivatives with respect to mass
    pub dr_dm: f64,
    pub dp_dm: f64,
    pub dl_dm: f64,
    pub dt_dm: f64,
}

#[derive(Debug, Clone)]
pub struct Star {
    pub shells: Vec<Shell>,

    // Initial / surface mass fractions
    pub x: f64, // Hydrogren
    pub y: f64, // Helium
    pub z: f64, // Metallicity

    // Per-shell hydrogen mass fraction on a fixed Lagrangian mass grid (COMP_GRID points,
    // index i corresponds to enclosed mass i * mass / COMP_GRID). This is what changes as
    // the star ages: the core burns first.
    pub comp_x: Vec<f64>,

    // Converged solution parameters (the four unknowns of the fitting method), reused as
    // the warm-start guess for the next solve.
    pub p_c: f64, // central pressure
    pub t_c: f64, // central temperature
    pub r_s: f64, // surface radius
    pub l_s: f64, // surface luminosity

    // Properties
    pub age: f64,  // Age in years
    pub mass: f64, // Mass in grams
}

impl Shell {
    // Update derived properties for this shell using the local hydrogen fraction `x`
    // (helium taken as 1 - x - z to conserve mass). The mass derivatives are recomputed by
    // the integrator; this fills in rho/kappa/epsilon/nabla for inspection and burning.
    pub fn update_derivatives(&mut self, x: f64, z: f64) {
        let (dr, dp, dl, dt) = get_derivatives(self.m, self.r, self.p, self.l, self.t, x, z);
        self.dr_dm = dr;
        self.dp_dm = dp;
        self.dl_dm = dl;
        self.dt_dm = dt;

        let y = (1.0 - x - z).max(0.0);
        let mu = mean_molecular_weight(x, y, z);
        self.rho = calculate_density(self.p, self.t, mu);
        self.kappa = calculate_opacity(self.rho, self.t, x, z);
        self.epsilon = calculate_epsilon(self.rho, self.t, x);
        // Actual temperature gradient nabla = dlnT/dlnP, recovered from the derivatives.
        self.nabla = if dp != 0.0 && self.t != 0.0 {
            (dt * self.p) / (dp * self.t)
        } else {
            0.0
        };
        self.x = x;
    }
}

impl Star {
    pub fn new(mass_solar: f64, x: f64, y: f64, z: f64) -> Self {
        // The star starts chemically homogeneous: every grid point has the surface composition.
        let comp_x = vec![x; COMP_GRID + 1];
        Star {
            shells: Vec::new(),
            x,
            y,
            z,
            comp_x,
            // Solar-like initial guesses (CGS)
            p_c: 2.4e17,
            t_c: 1.5e7,
            r_s: 6.957e10,
            l_s: 3.828e33,
            age: 0.0,
            mass: mass_solar * M_SUN,
        }
    }

    // Hydrogen fraction at enclosed mass `m`, read off the fixed composition grid.
    fn x_at_mass(&self, m: f64) -> f64 {
        let dm0 = self.mass / (COMP_GRID as f64);
        let idx = (m / dm0).round() as usize;
        self.comp_x[idx.min(self.comp_x.len() - 1)]
    }

    // Initialise the central shell with the standard power-series expansion about m = 0,
    // which removes the coordinate singularity at the centre.
    fn build_core(&self, p_c: f64, t_c: f64, dm: f64) -> Shell {
        let x0 = self.x_at_mass(0.0);
        let y0 = (1.0 - x0 - self.z).max(0.0);
        let mu = mean_molecular_weight(x0, y0, self.z);
        let rho_c = calculate_density(p_c, t_c, mu);
        let eps_c = calculate_epsilon(rho_c, t_c, x0);

        let r_0 = ((3.0 * dm) / (4.0 * PI * rho_c)).powf(1.0 / 3.0);
        let p_0 = p_c
            - (3.0 * G / (8.0 * PI))
                * ((4.0 * PI * rho_c / 3.0).powf(4.0 / 3.0))
                * dm.powf(2.0 / 3.0);
        let l_0 = eps_c * dm;

        let mut s = Shell {
            m: dm,
            r: r_0,
            p: p_0,
            l: l_0,
            t: t_c,
            x: x0,
            ..Default::default()
        };
        s.update_derivatives(x0, self.z);
        s
    }

    // Build the outer (surface) shell from trial values of the surface radius and luminosity.
    // The surface temperature follows from the black-body / Eddington photosphere,
    //   L = 4 pi R^2 sigma T_eff^4,
    // and the surface pressure from the optical-depth-2/3 photospheric condition,
    //   P_s = (2/3) g / kappa,
    // using electron-scattering opacity (a stable, well-behaved choice at the photosphere).
    fn build_surface(&self, r_s: f64, l_s: f64) -> Shell {
        let x = self.x_at_mass(self.mass);
        let t_eff = (l_s / (4.0 * PI * r_s.powi(2) * SIGMA)).powf(0.25);
        let g = G * self.mass / r_s.powi(2);
        let kappa_es = 0.2 * (1.0 + x);
        let p_s = (2.0 / 3.0) * g / kappa_es;

        let mut s = Shell {
            m: self.mass,
            r: r_s,
            p: p_s,
            l: l_s,
            t: t_eff,
            x,
            ..Default::default()
        };
        s.update_derivatives(x, self.z);
        s
    }

    // Adaptively integrate the structure equations from `start` to enclosed mass `m_target`
    // (in either direction). The step size is chosen so the fractional change in r, P and T
    // stays below STEP_TOL, which is what lets us march stably through both the smooth deep
    // interior and the stiff, radius-extended surface layers. Optionally records every shell.
    fn integrate_adaptive(
        &self,
        start: Shell,
        m_target: f64,
        mut collect: Option<&mut Vec<Shell>>,
    ) -> Shell {
        let mut s = start;
        if let Some(c) = collect.as_deref_mut() {
            c.push(s);
        }

        let dir = (m_target - s.m).signum();
        let dm_max = 0.05 * self.mass;

        let mut guard = 0;
        while (m_target - s.m) * dir > 0.0 && guard < 2_000_000 {
            guard += 1;

            let x = self.x_at_mass(s.m);
            let (dr, dp, _dl, dt) = get_derivatives(s.m, s.r, s.p, s.l, s.t, x, self.z);

            // Step from the fastest-changing of r, P, T so no variable moves by more than
            // STEP_TOL in a step. This must be allowed to become very small in the stiff
            // surface layers (which carry a negligible mass fraction) and grows naturally
            // into the smooth deep interior.
            let mut dm = dm_max;
            if dr.abs() > 0.0 {
                dm = dm.min(STEP_TOL * s.r.abs() / dr.abs());
            }
            if dp.abs() > 0.0 {
                dm = dm.min(STEP_TOL * s.p.abs() / dp.abs());
            }
            if dt.abs() > 0.0 {
                dm = dm.min(STEP_TOL * s.t.abs() / dt.abs());
            }
            dm = dm.min((m_target - s.m).abs());

            s = rk4_step(&s, dir * dm, x, self.z);

            if let Some(c) = collect.as_deref_mut() {
                c.push(s);
            }

            // Bail out of any numerical blow-up rather than spreading NaNs into the solver.
            if !s.p.is_finite() || !s.t.is_finite() || !s.r.is_finite() || s.p <= 0.0 || s.t <= 0.0
            {
                break;
            }
        }

        s
    }

    // Integrate outward from the centre to the fitting mass.
    fn integrate_out(&self, p_c: f64, t_c: f64, m_fit: f64, collect: Option<&mut Vec<Shell>>) -> Shell {
        let m_core = 1.0e-8 * self.mass;
        let core = self.build_core(p_c, t_c, m_core);
        self.integrate_adaptive(core, m_fit, collect)
    }

    // Integrate inward from the surface to the fitting mass.
    fn integrate_in(&self, r_s: f64, l_s: f64, m_fit: f64, collect: Option<&mut Vec<Shell>>) -> Shell {
        let surface = self.build_surface(r_s, l_s);
        self.integrate_adaptive(surface, m_fit, collect)
    }

    // The four fitting residuals: outward and inward integrations must agree on r, P, L and T
    // at the fitting point. Each is normalised by a characteristic scale so the Newton system
    // is well conditioned. A correct (p_c, t_c, r_s, l_s) drives all four to zero.
    fn fitting_residuals(&self, u: [f64; 4], m_fit: f64) -> [f64; 4] {
        let [p_c, t_c, r_s, l_s] = u;
        let out = self.integrate_out(p_c, t_c, m_fit, None);
        let inn = self.integrate_in(r_s, l_s, m_fit, None);
        [
            (out.r - inn.r) / r_s,
            (out.p - inn.p) / p_c,
            (out.l - inn.l) / l_s,
            (out.t - inn.t) / t_c,
        ]
    }

    // Solve the stellar structure as a two-point boundary value problem using the fitting
    // method: integrate outward from the core and inward from the surface (each in its
    // numerically stable direction) and require the two solutions to match at a fitting
    // point. A 4-D Newton-Raphson iteration adjusts (p_c, t_c, r_s, l_s) until they do.
    //
    // This replaces the original single un-iterated outward sweep, which could neither
    // satisfy the surface boundary conditions nor cope with the exponential error growth
    // that makes pure outward shooting unstable.
    pub fn solve_structure(&mut self) {
        let m_fit = 0.5 * self.mass;
        let mut u = [self.p_c, self.t_c, self.r_s, self.l_s];

        let norm = |f: &[f64; 4]| f.iter().map(|v| v * v).sum::<f64>().sqrt();
        let mut f = self.fitting_residuals(u, m_fit);

        for _ in 0..60 {
            if norm(&f) < 1e-6 {
                break;
            }

            // Finite-difference Jacobian J[i][j] = d f_i / d u_j
            let mut j = [[0.0f64; 4]; 4];
            for k in 0..4 {
                let h = u[k].abs() * 1e-4;
                let mut up = u;
                up[k] += h;
                let fp = self.fitting_residuals(up, m_fit);
                for i in 0..4 {
                    j[i][k] = (fp[i] - f[i]) / h;
                }
            }

            // Solve J * delta = -f
            let neg_f = [-f[0], -f[1], -f[2], -f[3]];
            let delta = match solve4(j, neg_f) {
                Some(d) => d,
                None => break, // Singular Jacobian
            };

            // Damped line search: shrink the step until the residual actually decreases and
            // the parameters stay physical (positive).
            let mut lambda = 1.0;
            let mut improved = false;
            for _ in 0..40 {
                let mut u_new = u;
                for i in 0..4 {
                    u_new[i] += lambda * delta[i];
                }
                if u_new.iter().all(|v| *v > 0.0) {
                    let f_new = self.fitting_residuals(u_new, m_fit);
                    if norm(&f_new) < norm(&f) {
                        u = u_new;
                        f = f_new;
                        improved = true;
                        break;
                    }
                }
                lambda *= 0.5;
            }

            if !improved {
                break;
            }
        }

        self.p_c = u[0];
        self.t_c = u[1];
        self.r_s = u[2];
        self.l_s = u[3];

        // Rebuild the full structure profile (centre -> fit, then surface -> fit reversed).
        let mut out_shells = Vec::new();
        self.integrate_out(u[0], u[1], m_fit, Some(&mut out_shells));
        let mut in_shells = Vec::new();
        self.integrate_in(u[2], u[3], m_fit, Some(&mut in_shells));
        in_shells.reverse();
        out_shells.extend(in_shells);
        self.shells = out_shells;
    }

    // Energy generation rate at enclosed mass `m`, interpolated from the solved profile.
    fn epsilon_at_mass(&self, m: f64) -> f64 {
        if self.shells.is_empty() {
            return 0.0;
        }
        // shells are ordered by increasing mass
        match self
            .shells
            .binary_search_by(|s| s.m.partial_cmp(&m).unwrap())
        {
            Ok(i) => self.shells[i].epsilon,
            Err(i) => {
                if i == 0 {
                    self.shells[0].epsilon
                } else if i >= self.shells.len() {
                    self.shells[self.shells.len() - 1].epsilon
                } else {
                    let a = &self.shells[i - 1];
                    let b = &self.shells[i];
                    let frac = (m - a.m) / (b.m - a.m).max(1.0);
                    a.epsilon + frac * (b.epsilon - a.epsilon)
                }
            }
        }
    }

    // Advance the chemical composition by one timestep and re-solve the structure.
    // Hydrogen is consumed where energy is generated: dX/dt = -epsilon / Q_H. As the core
    // hydrogen is depleted the mean molecular weight rises, the core contracts and heats and
    // the luminosity climbs -- the star evolves along the main sequence. This replaces the
    // original `evolve_timestep`, which only incremented the age.
    pub fn evolve_timestep(&mut self, dt_years: f64) {
        let dt_sec = dt_years * SEC_PER_YEAR;
        let dm0 = self.mass / (COMP_GRID as f64);

        for i in 0..self.comp_x.len() {
            let m = (i as f64) * dm0;
            let eps = self.epsilon_at_mass(m);
            let dx = eps * dt_sec / Q_H;
            self.comp_x[i] = (self.comp_x[i] - dx).max(0.0);
        }

        self.age += dt_years;
        self.solve_structure();
    }

    pub fn luminosity(&self) -> f64 {
        self.l_s
    }

    pub fn radius(&self) -> f64 {
        self.r_s
    }
}

// Solve a 4x4 linear system A x = b by Gaussian elimination with partial pivoting.
fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        // Partial pivot
        let mut piv = col;
        for row in (col + 1)..4 {
            if a[row][col].abs() > a[piv][col].abs() {
                piv = row;
            }
        }
        if a[piv][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);

        // Eliminate below
        for row in (col + 1)..4 {
            let factor = a[row][col] / a[col][col];
            for c in col..4 {
                a[row][c] -= factor * a[col][c];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution
    let mut x = [0.0f64; 4];
    for row in (0..4).rev() {
        let mut sum = b[row];
        for c in (row + 1)..4 {
            sum -= a[row][c] * x[c];
        }
        x[row] = sum / a[row][row];
    }
    Some(x)
}
