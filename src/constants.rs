/*
Important constants, all in CGS units
*/

pub const G: f64 = 6.67430e-8; // Gravitational constant
pub const A: f64 = 7.5657e-15; // Radiation density constant
pub const C: f64 = 2.99792e10; // Speed of light
pub const K_B: f64 = 1.380649e-16; // Boltzmann constant
pub const M_H: f64 = 1.6735e-24; // Hydrogen atom mass
pub const M_SUN: f64 = 1.989e33; // Mass of the sun
pub const SIGMA: f64 = 5.6704e-5; // Stefan-Boltzmann constant (= a*c/4)

// Nuclear energy released per gram of hydrogen fused into helium.
// ~0.7% of the rest-mass energy: 0.007 * c^2 ~ 6.3e18 erg/g
pub const Q_H: f64 = 6.3e18;

// Seconds in a year, for converting evolutionary timesteps into burn rates
pub const SEC_PER_YEAR: f64 = 3.156e7;
