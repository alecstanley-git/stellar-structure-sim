mod constants;
mod parameters;
mod physics;
mod star;

use physics::{calculate_density, mean_molecular_weight};
use star::Star;

use clap::Parser;
use serde::Serialize;

const R_SUN: f64 = 6.957e10; // cm
const L_SUN: f64 = 3.828e33; // erg/s

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(long, default_value_t = 1.0)]
    mass: f64,
    #[arg(long, default_value_t = 0.70)]
    x: f64,
    #[arg(long, default_value_t = 0.28)]
    y: f64,
    #[arg(long, default_value_t = 0.02)]
    z: f64,
    #[arg(long, default_value_t = 10.0)]
    max_age: f64,
    #[arg(long)]
    json: bool,
}

#[derive(Serialize)]
struct StateRecord {
    age_gyr: f64,
    l_lsun: f64,
    r_rsun: f64,
    teff_k: f64,
    pc_cgs: f64,
    x_core: f64,
}

fn record_state(star: &Star) -> StateRecord {
    let core = &star.shells[0];
    let surface = &star.shells[star.shells.len() - 1];
    let teff = (surface.l / (4.0 * std::f64::consts::PI * surface.r.powi(2) * 5.67e-5)).powf(0.25);

    StateRecord {
        age_gyr: star.age / 1e9,
        l_lsun: surface.l / L_SUN,
        r_rsun: surface.r / R_SUN,
        teff_k: teff,
        pc_cgs: core.p,
        x_core: core.x,
    }
}

fn main() {
    let args = Args::parse();
    let mut star = Star::new(args.mass, args.x, args.y, args.z);
    let mut history = Vec::new();

    if !args.json {
        println!(
            "Composition: X = {:.2}, Y = {:.2}, Z = {:.2}",
            star.x, star.y, star.z
        );
        println!("--- Zero-Age Main Sequence ---");
    }

    star.solve_structure(0.0);

    if args.json {
        history.push(record_state(&star));
    } else {
        print_state(&star);
        println!("--- Main-Sequence Evolution ---");
    }

    let mut timestep = 0.5e9 * args.mass.powf(-2.5); // years

    while star.age / 1e9 < args.max_age {
        // Prevent overshooting max_age
        if (star.age + timestep) / 1e9 > args.max_age {
            timestep = args.max_age * 1e9 - star.age;
        }

        star.evolve_timestep(timestep);
        if args.json {
            history.push(record_state(&star));
        } else {
            print_state(&star);
        }
    }

    if args.json {
        let json_output = serde_json::to_string(&history).unwrap();
        println!("{}", json_output);
    }
}

pub fn print_state(star: &Star) {
    let s0 = &star.shells[0];
    let y0 = (1.0 - s0.x - star.z).max(0.0);
    let mu0 = mean_molecular_weight(s0.x, y0, star.z);
    let rho_c = calculate_density(s0.p, s0.t, mu0);
    println!(
        "age = {:5.2} Gyr | L = {:.3} L_sun | R = {:.3} R_sun | P_c = {:.3e} | T_c = {:.3e} K | X_core = {:.3} | rho_c = {:.1} | r_core = {:.3e}",
        star.age / 1e9,
        star.luminosity() / L_SUN,
        star.radius() / R_SUN,
        s0.p,
        s0.t,
        s0.x,
        rho_c,
        s0.r
    );

    if star.age < 0.1e9 || star.age > 9.9e9 {
        let n = star.shells.len();
        println!(
            "  Shell 500: P={:.3e}, R={:.3e}, T={:.3e}, M={:.3e}",
            star.shells[500].p, star.shells[500].r, star.shells[500].t, star.shells[500].m
        );
        println!(
            "  Shell 1500: P={:.3e}, R={:.3e}, T={:.3e}, M={:.3e}",
            star.shells[1500].p, star.shells[1500].r, star.shells[1500].t, star.shells[1500].m
        );
    }
}
