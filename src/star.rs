use std::f64::consts::PI;

use crate::{
    constants::{G, M_SUN, Q_H, SEC_PER_YEAR, SIGMA},
    parameters::COMP_GRID,
    physics::{
        calculate_density, calculate_epsilon, calculate_epsilon_grav, calculate_nabla,
        calculate_opacity, mean_molecular_weight,
    },
};

const R_SCALE: f64 = 6.957e10;
const P_SCALE: f64 = 1e17;
const L_SCALE: f64 = 3.828e33;
const T_SCALE: f64 = 1e7;

#[derive(Debug, Clone, Copy, Default)]
pub struct Shell {
    pub m: f64,
    pub r: f64,
    pub p: f64,
    pub l: f64,
    pub t: f64,
    pub x: f64,
    pub y: f64,

    pub p_old: f64,
    pub t_old: f64,
}

#[derive(Debug, Clone)]
pub struct Star {
    pub shells: Vec<Shell>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub age: f64,
    pub mass: f64,
}

impl Star {
    pub fn new(mass_solar: f64, x: f64, y: f64, z: f64) -> Self {
        let n = COMP_GRID;
        let mass = mass_solar * M_SUN;
        let mut shells = Vec::with_capacity(n);

        let min_m = 1.0e-8 * mass;
        for i in 0..n {
            let f = (i as f64) / ((n - 1) as f64);
            let m_frac = if f < 0.5 {
                0.5 * (2.0 * f).powf(2.5)
            } else {
                1.0 - 0.5 * (2.0 * (1.0 - f)).powf(6.0)
            };
            let m = min_m + m_frac * (mass - min_m);

            // Initial guess using homology relations
            let frac = m / mass;
            let r_scale_factor = mass_solar.powf(0.75);
            let p_scale_factor = mass_solar.powf(-0.8);
            let t_scale_factor = mass_solar.powf(0.3);
            let l_scale_factor = mass_solar.powf(3.5);

            let r = 0.9 * 6.957e10 * r_scale_factor * frac.powf(1.0 / 3.0);
            let p = 2.4e17 * p_scale_factor * (1.0 - frac.powf(2.0 / 3.0)).max(1e-10) + 1e4;
            let t = 1.5e7 * t_scale_factor * (1.0 - frac.powf(2.0 / 3.0)).max(1e-10)
                + 5770.0 * mass_solar.powf(0.1);
            let l = 0.9 * 3.828e33 * l_scale_factor * frac;

            shells.push(Shell {
                m,
                r,
                p,
                l,
                t,
                x,
                y,
                p_old: p,
                t_old: t,
            });
        }

        let mut star = Star {
            shells,
            x,
            y,
            z,
            age: 0.0,
            mass,
        };

        star.solve_structure(0.0);

        for i in 0..n {
            star.shells[i].p_old = star.shells[i].p;
            star.shells[i].t_old = star.shells[i].t;
        }

        star
    }

    fn calculate_residuals(&self, dt: f64) -> Vec<f64> {
        let n = self.shells.len();
        let mut f = vec![0.0; 4 * n];
        for i in 0..n {
            self.calc_eqs(dt, &mut f, i);
        }
        f
    }

    fn calc_eqs(&self, dt: f64, f: &mut [f64], i: usize) {
        let n = self.shells.len();

        // Inner boundary (eq 0, 1)
        if i == 0 {
            let s0 = &self.shells[0];
            let mu0 = mean_molecular_weight(s0.x, s0.y, self.z);
            let rho0 = calculate_density(s0.p, s0.t, mu0).max(1e-10);
            let (eps_h, eps_he) = calculate_epsilon(rho0, s0.t, s0.x, s0.y, self.z);
            let eps_nuc0 = eps_h + eps_he;
            let eps_grav0 = calculate_epsilon_grav(s0.t, s0.t_old, s0.p, s0.p_old, mu0, dt);

            f[0] = (s0.r.powi(3) - (3.0 * s0.m) / (4.0 * PI * rho0)) / s0.r.powi(3).max(1e-20);
            f[1] = (s0.l - (eps_nuc0 + eps_grav0) * s0.m) / (s0.l.abs() + 1e-4 * L_SCALE);
        }

        // Structure equations (eq 4*j+2 .. 4*j+5)
        // We evaluate for j = i-1 and j = i
        for j in i.saturating_sub(1)..=i {
            if j >= n - 1 {
                continue;
            }
            let s1 = &self.shells[j];
            let s2 = &self.shells[j + 1];

            let dm = s2.m - s1.m;
            let r_bar = 0.5 * (s1.r + s2.r).max(1e-10);
            let p_bar = 0.5 * (s1.p + s2.p).max(1e-10);
            let l_bar = 0.5 * (s1.l + s2.l);
            let t_bar = 0.5 * (s1.t + s2.t).max(1e-10);
            let x_bar = 0.5 * (s1.x + s2.x);
            let t_old_bar = 0.5 * (s1.t_old + s2.t_old).max(1e-10);
            let p_old_bar = 0.5 * (s1.p_old + s2.p_old).max(1e-10);
            let m_bar = 0.5 * (s1.m + s2.m);

            let y_bar = 0.5 * (s1.y + s2.y);
            let mu = mean_molecular_weight(x_bar, y_bar, self.z);
            let rho = calculate_density(p_bar, t_bar, mu).max(1e-10);
            let kappa = calculate_opacity(rho, t_bar, x_bar, self.z);
            let (eps_h, eps_he) = calculate_epsilon(rho, t_bar, x_bar, y_bar, self.z);
            let eps_nuc = eps_h + eps_he;
            let eps_grav = calculate_epsilon_grav(t_bar, t_old_bar, p_bar, p_old_bar, mu, dt);

            // Exact integral for r^3 = r^3 + 3dm / (4 pi rho)
            let f_r3 = 3.0 / (4.0 * PI * rho);
            // Analytic integration for P using m^{2/3} for inner shells
            let f_p = if j < 5 {
                -(3.0 * G / (8.0 * PI))
                    * (4.0 * PI * rho / 3.0).powf(4.0 / 3.0)
                    * (s2.m.powf(2.0 / 3.0) - s1.m.powf(2.0 / 3.0))
                    / dm
            } else {
                -(G * m_bar) / (4.0 * PI * r_bar.powi(4))
            };
            let f_l = eps_nuc + eps_grav;

            // Mixing Length Theory parameter
            let alpha_mlt = 1.5;
            let nabla =
                calculate_nabla(t_bar, p_bar, m_bar, r_bar, l_bar, rho, kappa, mu, alpha_mlt);
            let f_t = -(G * m_bar * t_bar * nabla) / (4.0 * PI * r_bar.powi(4) * p_bar);

            let eq_idx = 4 * j + 2;
            f[eq_idx] = (s2.r.powi(3) - s1.r.powi(3) - dm * f_r3) / r_bar.powi(3);
            f[eq_idx + 1] = (s2.p - s1.p - dm * f_p) / p_bar;
            f[eq_idx + 2] = (s2.l - s1.l - dm * f_l) / (l_bar.abs() + 1e-4 * L_SCALE);
            f[eq_idx + 3] = (s2.t - s1.t - dm * f_t) / t_bar;
        }

        // Outer boundary
        if i == n - 1 {
            let s_n = &self.shells[n - 1];
            let mu_n = mean_molecular_weight(s_n.x, s_n.y, self.z);
            let rho_n = calculate_density(s_n.p, s_n.t, mu_n).max(1e-10);
            let kappa_n = calculate_opacity(rho_n, s_n.t, s_n.x, self.z);

            f[4 * n - 2] = (s_n.l - 4.0 * PI * s_n.r.powi(2) * SIGMA * s_n.t.powi(4))
                / (s_n.l.abs() + 1e-4 * L_SCALE);
            f[4 * n - 1] = (s_n.p - (2.0 / 3.0) * (G * self.mass) / (s_n.r.powi(2) * kappa_n))
                / s_n.p.max(1e-10);
        }
    }

    pub fn solve_structure(&mut self, dt: f64) {
        let n = self.shells.len();
        let mut converged = false;
        let norm = |f: &[f64]| (f.iter().map(|v| v * v).sum::<f64>() / (f.len() as f64)).sqrt();
        let mut last_norm = std::f64::MAX;

        for iter in 0..200 {
            let f = self.calculate_residuals(dt);
            let current_norm = norm(&f);

            if current_norm < 1e-4 {
                converged = true;
                break;
            }
            if iter > 10 && (last_norm - current_norm).abs() / current_norm < 1e-4 {
                // If it stops improving, accept it
                converged = true;
                break;
            }
            last_norm = current_norm;

            // Construct block tridiagonal Jacobian
            let mut a_blocks = vec![[[0.0; 4]; 4]; n];
            let mut b_blocks = vec![[[0.0; 4]; 4]; n];
            let mut c_blocks = vec![[[0.0; 4]; 4]; n];
            let mut d_vecs = vec![[0.0; 4]; n];

            for i in 0..n {
                for eq in 0..4 {
                    d_vecs[i][eq] = -f[4 * i + eq];
                }
            }

            // Perturb each variable to compute Jacobian
            for i in 0..n {
                let vars = [
                    self.shells[i].r,
                    self.shells[i].p,
                    self.shells[i].l,
                    self.shells[i].t,
                ];
                let scales = [R_SCALE, P_SCALE, L_SCALE, T_SCALE];

                for var_idx in 0..4 {
                    let h = (vars[var_idx].abs() * 1e-6).max(scales[var_idx] * 1e-8);

                    // Modify var
                    let orig = vars[var_idx];
                    match var_idx {
                        0 => self.shells[i].r += h,
                        1 => self.shells[i].p += h,
                        2 => self.shells[i].l += h,
                        3 => self.shells[i].t += h,
                        _ => {}
                    }

                    let mut f_new = f.clone();
                    self.calc_eqs(dt, &mut f_new, i);

                    // Restore var
                    match var_idx {
                        0 => self.shells[i].r = orig,
                        1 => self.shells[i].p = orig,
                        2 => self.shells[i].l = orig,
                        3 => self.shells[i].t = orig,
                        _ => {}
                    }

                    // Fill B_i
                    for eq in 0..4 {
                        b_blocks[i][eq][var_idx] = (f_new[4 * i + eq] - f[4 * i + eq]) / h;
                    }

                    // Fill A_i
                    if i < n - 1 {
                        for eq in 0..4 {
                            a_blocks[i + 1][eq][var_idx] =
                                (f_new[4 * (i + 1) + eq] - f[4 * (i + 1) + eq]) / h;
                        }
                    }

                    // Fill C_i
                    if i > 0 {
                        for eq in 0..4 {
                            c_blocks[i - 1][eq][var_idx] =
                                (f_new[4 * (i - 1) + eq] - f[4 * (i - 1) + eq]) / h;
                        }
                    }
                }
            }

            // Solve block tridiagonal system
            if let Some(delta) = solve_block_tridiagonal(&a_blocks, &b_blocks, &c_blocks, &d_vecs) {
                let mut max_frac = 0.0_f64;
                for i in 0..n {
                    let s = &self.shells[i];
                    max_frac = max_frac
                        .max(delta[i][0].abs() / s.r.max(R_SCALE * 1e-4))
                        .max(delta[i][1].abs() / s.p.max(P_SCALE * 1e-4))
                        .max(delta[i][2].abs() / (s.l.abs() + L_SCALE * 1e-4))
                        .max(delta[i][3].abs() / s.t.max(T_SCALE * 1e-4));
                }

                let lambda = if max_frac > 0.3 { 0.3 / max_frac } else { 1.0 };

                for i in 0..n {
                    self.shells[i].r += lambda * delta[i][0];
                    self.shells[i].p += lambda * delta[i][1];
                    self.shells[i].l += lambda * delta[i][2];
                    self.shells[i].t += lambda * delta[i][3];
                    // Constrain to positive values physically
                    self.shells[i].r = self.shells[i].r.max(1e2);
                    self.shells[i].p = self.shells[i].p.max(1e2);
                    self.shells[i].t = self.shells[i].t.max(10.0);
                }
            } else {
                break;
            }
        }

        if !converged {
            println!("Warning: Solver failed to converge!");
        }
    }

    pub fn evolve_timestep(&mut self, dt_years: f64) {
        let dt_sec = dt_years * SEC_PER_YEAR;
        let q_he = 6.0e17; // Energy per gram from Triple-Alpha roughly

        let mut convective_flags = vec![false; self.shells.len()];

        // Burn hydrogen and helium explicitly, and check convective regions
        for i in 0..self.shells.len() {
            let (p, t, x, y, m, r, l) = {
                let s = &self.shells[i];
                (s.p, s.t, s.x, s.y, s.m, s.r, s.l)
            };

            let mu = mean_molecular_weight(x, y, self.z);
            let rho = calculate_density(p, t, mu).max(1e-10);
            let (eps_h, eps_he) = calculate_epsilon(rho, t, x, y, self.z);

            let dx = eps_h * dt_sec / Q_H;
            let dy = eps_he * dt_sec / q_he;

            let s = &mut self.shells[i];

            let real_dx = dx.min(s.x);
            s.x -= real_dx;
            s.y += real_dx; // H burning creates He

            let real_dy = dy.min(s.y);
            s.y -= real_dy;

            s.p_old = p;
            s.t_old = t;

            // Check convection for mixing
            let kappa = calculate_opacity(rho, t, x, self.z);
            let nabla_rad = (3.0 * kappa * l * p)
                / (16.0 * PI * 7.56e-15 * 3e10 * crate::constants::G * m * t.powi(4));
            let nabla_ad = 0.4;
            convective_flags[i] = nabla_rad > nabla_ad;
        }

        // Convective Mixing
        let mut start_idx = 0;
        while start_idx < self.shells.len() {
            if convective_flags[start_idx] {
                let mut end_idx = start_idx;
                while end_idx < self.shells.len() && convective_flags[end_idx] {
                    end_idx += 1;
                }

                // Mix composition between start_idx and end_idx
                let mut total_mass = 0.0;
                let mut total_x_mass = 0.0;
                let mut total_y_mass = 0.0;

                for i in start_idx..end_idx {
                    let dm = if i == 0 {
                        self.shells[0].m
                    } else {
                        self.shells[i].m - self.shells[i - 1].m
                    };
                    total_mass += dm;
                    total_x_mass += self.shells[i].x * dm;
                    total_y_mass += self.shells[i].y * dm;
                }

                let avg_x = total_x_mass / total_mass.max(1e-30);
                let avg_y = total_y_mass / total_mass.max(1e-30);

                for i in start_idx..end_idx {
                    self.shells[i].x = avg_x;
                    self.shells[i].y = avg_y;
                }

                start_idx = end_idx;
            } else {
                start_idx += 1;
            }
        }

        self.age += dt_years;
        self.solve_structure(dt_sec);
    }

    pub fn luminosity(&self) -> f64 {
        self.shells.last().unwrap().l
    }

    pub fn radius(&self) -> f64 {
        self.shells.last().unwrap().r
    }
}

// Block tridiagonal solver
fn solve_block_tridiagonal(
    a: &[[[f64; 4]; 4]],
    b: &[[[f64; 4]; 4]],
    c: &[[[f64; 4]; 4]],
    d: &[[f64; 4]],
) -> Option<Vec<[f64; 4]>> {
    let n = b.len();
    let mut c_prime = vec![[[0.0; 4]; 4]; n];
    let mut d_prime = vec![[0.0; 4]; n];

    let inv_b0 = invert_4x4(b[0])?;
    c_prime[0] = mat_mul_4x4(inv_b0, c[0]);
    d_prime[0] = mat_vec_mul_4x4(inv_b0, d[0]);

    for i in 1..n {
        let mut temp = b[i];
        let ac = mat_mul_4x4(a[i], c_prime[i - 1]);
        for r in 0..4 {
            for col in 0..4 {
                temp[r][col] -= ac[r][col];
            }
        }

        let inv_temp = invert_4x4(temp)?;
        if i < n - 1 {
            c_prime[i] = mat_mul_4x4(inv_temp, c[i]);
        }

        let mut vec_temp = d[i];
        let ad = mat_vec_mul_4x4(a[i], d_prime[i - 1]);
        for r in 0..4 {
            vec_temp[r] -= ad[r];
        }
        d_prime[i] = mat_vec_mul_4x4(inv_temp, vec_temp);
    }

    let mut x = vec![[0.0; 4]; n];
    x[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        let cx = mat_vec_mul_4x4(c_prime[i], x[i + 1]);
        for r in 0..4 {
            x[i][r] = d_prime[i][r] - cx[r];
        }
    }

    Some(x)
}

fn invert_4x4(m: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut inv = [[0.0; 4]; 4];
    for col in 0..4 {
        let mut b = [0.0; 4];
        b[col] = 1.0;
        if let Some(x) = solve4(m, b) {
            for row in 0..4 {
                inv[row][col] = x[row];
            }
        } else {
            return None;
        }
    }
    Some(inv)
}

fn mat_mul_4x4(m1: [[f64; 4]; 4], m2: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut res = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    res
}

fn mat_vec_mul_4x4(m: [[f64; 4]; 4], v: [f64; 4]) -> [f64; 4] {
    let mut res = [0.0; 4];
    for i in 0..4 {
        for j in 0..4 {
            res[i] += m[i][j] * v[j];
        }
    }
    res
}

fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
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

        for row in (col + 1)..4 {
            let factor = a[row][col] / a[col][col];
            let (left, right) = a.split_at_mut(row);
            for (target_val, &pivot_val) in right[0].iter_mut().zip(left[col].iter()).skip(col) {
                *target_val -= factor * pivot_val;
            }
            b[row] -= factor * b[col];
        }
    }

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
