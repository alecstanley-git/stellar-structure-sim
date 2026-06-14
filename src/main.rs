mod constants;
mod integrator;
mod parameters;
mod physics;
mod star;

use star::Star;

// Handy solar reference values for sanity-checking the output (CGS)
const R_SUN: f64 = 6.957e10; // cm
const L_SUN: f64 = 3.828e33; // erg/s

fn main() {
    let mut star = Star::new(1.0, 0.70, 0.28, 0.02);
    println!(
        "Composition: X = {:.2}, Y = {:.2}, Z = {:.2}",
        star.x, star.y, star.z
    );

    // Solve the zero-age structure as a two-point boundary value problem (fitting method).
    // Initial guesses for (P_c, T_c, R_s, L_s) are seeded in Star::new from solar values.
    star.solve_structure();

    println!("--- Zero-Age Main Sequence ---");
    print_state(&star);

    // Evolve along the main sequence, burning core hydrogen each step.
    println!("--- Main-Sequence Evolution ---");
    let timestep = 1.0e9; // years
    for _ in 0..10 {
        star.evolve_timestep(timestep);
        print_state(&star);
    }
}

fn print_state(star: &Star) {
    println!(
        "age = {:5.2} Gyr | L = {:.3} L_sun | R = {:.3} R_sun | P_c = {:.3e} | T_c = {:.3e} K | X_core = {:.3}",
        star.age / 1.0e9,
        star.luminosity() / L_SUN,
        star.radius() / R_SUN,
        star.p_c,
        star.t_c,
        star.comp_x[0],
    );
}
