/// Compute one coefficient for the power-law distribution, using normalized
/// parameters that more directly speak to the variables that control the
/// problem.

extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI, Coefficient, Stokes};

fn main() {
    const S: f64 = 5.476500e5;
    const THETA: f64 = 1.381850e0;
    const P: f64 = 2.394880e0;
    const COEFF: Coefficient = Coefficient::Emission(Stokes::I);
    const SYMPHONY_JI: f64 = 1.3415653497801385e-30;

    const NU: f64 = 1e9;
    const B: f64 = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * S);
    const N_E: f64 = 1.;

    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1000.;
    const GAMMA_CUTOFF: f64 = 1e7;

    let ji = rimphony::PowerLawDistribution::new(P)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .finish(COEFF, NU, B, N_E, THETA)
        .compute();

    println!("Symphony j_I: {:e}   Ours: {:e}", SYMPHONY_JI, ji);
}
