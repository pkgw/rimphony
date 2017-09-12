/// Compute one coefficient for the power-law distribution, using normalized
/// parameters that more directly speak to the variables that control the
/// problem.

extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI, Coefficient, Stokes};

fn main() {
    const S: f64 = 1.3136507771e2;
    const THETA: f64 = 8.7361125358e-1;
    const P: f64 = 1.5498441750e0;
    const COEFF: Coefficient = Coefficient::Absorption(Stokes::I);
    const SYMPHONY_VAL: f64 = 9.286451e-16;

    const NU: f64 = 1e9;
    const B: f64 = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * S);
    const N_E: f64 = 1.;

    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1000.;
    const GAMMA_CUTOFF: f64 = 1e7;

    let val = rimphony::PowerLawDistribution::new(P)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .finish(COEFF, NU, B, N_E, THETA)
        .compute();

    println!("Symphony: {:e}   Us: {:e}", SYMPHONY_VAL, val);
}
