/// Compute one coefficient for the power-law distribution, using the same
/// parameters that are fed directly into Symphony.

extern crate rimphony;

use rimphony::{Coefficient, Stokes, SynchrotronCalculator};

fn main() {
    const NU: f64 = 1e9;
    const B: f64 = 1e3;
    const N_E: f64 = 1.;
    const THETA: f64 = 0.9;
    const P: f64 = 2.5;
    const SYMPHONY_JI: f64 = 2.64399749412774e-21;

    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1e12;
    const GAMMA_CUTOFF: f64 = 1e10;

    let ji = rimphony::PowerLawDistribution::new(P)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .full_calculation()
        .compute_cgs(Coefficient::Emission, Stokes::I, NU, B, N_E, THETA);

    println!("Symphony j_I: {:e}   Ours: {:e}", SYMPHONY_JI, ji);
}
