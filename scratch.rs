use std::sync::Arc;
mod constants;
mod opacity;
mod parameters;
mod physics;
mod star;

fn main() {
    let table = opacity::OpacityTable::load_from_file("opacity_AGSS09_z0.02.txt").unwrap();
    let table_arc = Arc::new(table);
    let mut star = star::Star::new(1.0, 0.70, 0.28, 0.02, table_arc.clone());
    
    for i in 0..star.shells.len() {
        let s = &star.shells[i];
        let mu = physics::mean_molecular_weight(s.x, (1.0 - s.x - 0.02).max(0.0), 0.02);
        let rho = physics::calculate_density(s.p, s.t, mu);
        
        let k_new = physics::calculate_opacity(rho, s.t, s.x, 0.02, &table_arc);
        
        // Let's compute old_kappa
        let t6 = s.t / 1.0e6;
        let guillotine = 1.0 + (t6 * t6);
        let kappa_ff = 1.75e22 * (1.0 + s.x) * (1.0 - 0.02) * rho * s.t.powf(-3.5);
        let kappa_bf = 4.3e25 * 0.02 * (1.0 + s.x) * rho * s.t.powf(-3.5) / guillotine;
        let kappa_es = 0.2 * (1.0 + s.x);
        let kappa_interior = kappa_ff + kappa_bf + kappa_es;
        let kappa_h_minus = 2.5e-32 * (0.02 / 0.02) * rho.sqrt() * s.t.powi(9);
        let old_k = (kappa_h_minus * kappa_interior) / (kappa_h_minus + kappa_interior).max(1e-10);
        
        if i % 10 == 0 || i > star.shells.len() - 10 {
            println!("Shell {}: T={:.2e}, rho={:.2e}, old_K={:.4e}, new_K={:.4e}", i, s.t, rho, old_k, k_new);
        }
    }
}
