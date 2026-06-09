mod star;
mod vec3;

use star::Star;
use vec3::Vec3;

fn main() {
    let v1 = Vec3 {x: 0.0, y: 1.0, z: 0.0};
    let v2 = Vec3 {x: -1.0, y: 10.0, z: 5.0};
    // let v3 = v1 + v2;
    let v4 = v1 * v2;
    v4.print();
}
