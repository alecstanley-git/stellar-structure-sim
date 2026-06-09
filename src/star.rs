pub struct Star {
    age: f64,
    mass: f64
}

impl Star {
    pub fn new(age: f64, mass: f64) -> Self {
        Star { age, mass }
    }

    pub fn print(self) {
        println!("--- STAR PROPERTIES ---");
        println!("Age: {}", self.age.to_string());
        println!("Mass: {}", self.mass.to_string());
    }
}
