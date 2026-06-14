// Resolution of the fixed Lagrangian composition grid (hydrogen fraction vs. enclosed
// mass). The structure itself is integrated with an adaptive step size, but the chemical
// profile that evolves in time lives on this fixed mass grid.
pub const COMP_GRID: usize = 2000;

// Fractional change in (r, P, T) allowed per adaptive integration step.
pub const STEP_TOL: f64 = 0.01;
