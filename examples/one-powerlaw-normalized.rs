/// Compute one coefficient for the power-law distribution, using normalized
/// parameters that more directly speak to the variables that control the
/// problem.

extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI,
               Coefficient, Stokes, SynchrotronCalculator};

fn main() {
    const S: f64 = 1.0360583634e3;
    const THETA: f64 = 7.4017422303e-1;
    const P: f64 = 2.3306843452e0;
    const COEFF: Coefficient = Coefficient::Absorption;
    const STOKES: Stokes = Stokes::I;
    const SYMPHONY_VAL: f64 = 0.;

    const NU: f64 = 1e9;
    const B: f64 = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * S);
    const N_E: f64 = 1.;

    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1e12;
    const GAMMA_CUTOFF: f64 = 1e10;

    let val = rimphony::PowerLawDistribution::new(P)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .full_calculation()
        .compute_cgs(COEFF, STOKES, NU, B, N_E, THETA);

    let remove_units = match COEFF {
        Coefficient::Emission =>
            SPEED_LIGHT * THETA.cos().abs() / ((TWO_PI * ELECTRON_CHARGE).powi(2) * NU),
        Coefficient::Absorption =>
            -2. * MASS_ELECTRON * SPEED_LIGHT * NU * THETA.cos().abs() / (TWO_PI * ELECTRON_CHARGE).powi(2),
    };

    println!("Inner Symphony: {:e}   Us: {:e}", SYMPHONY_VAL * remove_units, val * remove_units);
    println!("Outer Symphony: {:e}   Us: {:e}", SYMPHONY_VAL, val);
}
