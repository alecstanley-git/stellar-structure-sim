use crate::parameters::INIT_SHELL_COUNT;

#[derive(Debug, Clone, Copy, Default)]
pub struct Shell {
    pub r: f64, // Radius
    pub p: f64, // Pressure
    pub l: f64, // Luminosity
    pub t: f64, // Temp
    pub dr_dm: f64,
    pub dp_dm: f64,
    pub dl_dm: f64,
    pub dt_dm: f64,
}

#[derive(Debug, Clone)]
pub struct Star {
    pub shells: Vec<Shell>,
}

impl Default for Star {
    fn default() -> Star {
        Star {
            shells: vec![Shell::default(); INIT_SHELL_COUNT],
        }
    }
}
