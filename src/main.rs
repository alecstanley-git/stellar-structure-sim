mod constants;
mod parameters;
mod physics;
mod star;
mod opacity;

use star::Star;
use opacity::OpacityTable;
use std::sync::Arc;

// Handy solar reference values for sanity-checking the output (CGS)
const R_SUN: f64 = 6.957e10; // cm
const L_SUN: f64 = 3.828e33; // erg/s

fn main() {
    println!("Loading opacity table...");
    let table = OpacityTable::load_from_file("opacity_AGSS09_z0.02.txt").expect("Failed to load opacity table");
    let table_arc = Arc::new(table);

    let mut star = Star::new(1.0, 0.70, 0.28, 0.02, table_arc);
    println!(
        "Composition: X = {:.2}, Y = {:.2}, Z = {:.2}",
        star.x, star.y, star.z
    );

    println!("--- Initial Guess (Old Physics) ---");
    star.aesopus_blend = 0.0;
    star.solve_structure(0.0);
    print_state(&star);

    println!("--- Relaxing into AESOPUS Opacities ---");
    for i in 1..=20 {
        star.aesopus_blend = (i as f64) / 20.0;
        star.solve_structure(0.0);
    }

    println!("--- Zero-Age Main Sequence (New Physics) ---");
    print_state(&star);

    // Evolve along the main sequence
    println!("--- Main-Sequence Evolution ---");
    let timestep = 0.5e9; // years
    for _ in 0..20 {
        star.evolve_timestep(timestep);
        print_state(&star);
    }
}

fn print_state(star: &Star) {
    let s_c = &star.shells[0];
    println!(
        "age = {:5.2} Gyr | L = {:.3} L_sun | R = {:.3} R_sun | P_c = {:.3e} | T_c = {:.3e} K | X_core = {:.3}",
        star.age / 1.0e9,
        star.luminosity() / L_SUN,
        star.radius() / R_SUN,
        s_c.p,
        s_c.t,
        s_c.x,
    );
}
