/// Compute one coefficient for the power-law distribution.
///
/// This uses the parameters from my first Symphony sample.

extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI, Coefficient, Stokes};

fn main() {
    const S: f64 = 5.39186e+07;
    const THETA: f64 = 1.29985;
    const P: f64 = 1.70877;
    const NU: f64 = 1e9;
    const N_E: f64 = 1.;
    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1000.;
    const GAMMA_CUTOFF: f64 = 1e7;
    const SYMPHONY_JI: f64 = 1.63868e-31;
    const B: f64 = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * S);

    let ji = rimphony::PowerLawDistribution::new(P)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .finish(Coefficient::Emission(Stokes::I), NU, B, N_E, THETA)
        .compute();

    println!("Symphony j_I: {}   Ours: {}", SYMPHONY_JI, ji);
}
