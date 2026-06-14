mod constants;
mod star;
mod parameters;

use star::Star;

fn main() {
    let star = Star::default();
    println!("{:#?}", star);
}
q