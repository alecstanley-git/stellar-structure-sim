mod constants;
mod parameters;
mod physics;
mod star;

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
    #[arg(long)]
    json: bool,
}

#[derive(Serialize)]
struct StateRecord {
    age_gyr: f64,
    l_lsun: f64,
    r_rsun: f64,
    tc_k: f64,
    pc_cgs: f64,
    x_core: f64,
}

fn record_state(star: &Star) -> StateRecord {
    let core = &star.shells[0];
    let surface = &star.shells[star.shells.len() - 1];
    StateRecord {
        age_gyr: star.age / 1e9,
        l_lsun: surface.l / L_SUN,
        r_rsun: surface.r / R_SUN,
        tc_k: core.t,
        pc_cgs: core.p,
        x_core: core.x,
    }
}

fn main() {
    let args = Args::parse();
    let mut star = Star::new(args.mass, args.x, args.y, args.z);
    let mut history = Vec::new();

    if !args.json {
        println!("Composition: X = {:.2}, Y = {:.2}, Z = {:.2}", star.x, star.y, star.z);
        println!("--- Zero-Age Main Sequence ---");
    }
    
    star.solve_structure(0.0);
    
    if args.json {
        history.push(record_state(&star));
    } else {
        print_state(&star);
        println!("--- Main-Sequence Evolution ---");
    }

    let timestep = 0.5e9; // years
    for _ in 0..20 {
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
