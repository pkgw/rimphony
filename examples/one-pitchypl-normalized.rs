// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Compute one coefficient for the pitchy power-law distribution, using
/// normalized parameters that more directly speak to the variables that
/// control the problem.

extern crate rimphony;

use rimphony::{Coefficient, Stokes, SynchrotronCalculator};

fn main() {
    const S: f64 = 8.0973407678629616e0;
    const THETA: f64 = 7.2687065355210786e-2;
    const P: f64 = 2.7273434060193211e0;
    const K: f64 = 2.7016346500930695e0;
    const COEFF: Coefficient = Coefficient::Faraday;
    const STOKES: Stokes = Stokes::Q;

    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1e12;
    const GAMMA_CUTOFF: f64 = 1e10;

    let val = rimphony::PitchyPowerLawDistribution::new(P, K)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .full_calculation()
        .compute_dimensionless(COEFF, STOKES, S, THETA);

    println!("{:.18e}", val);
}
