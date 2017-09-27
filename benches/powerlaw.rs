// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Crunch some numbers for in the powerlaw model.
///
/// Due to the way it’s designed, the benchmarker runs each "iteration" at
/// least 300 times, so this program can take a while to run. It adds up to
/// about 40 minutes on my machine.

#[macro_use] extern crate bencher;
extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI,
               Coefficient, Stokes, SynchrotronCalculator};
use bencher::Bencher;


// I designed a neat little Latin square system to try to run benchmarks
// across a moderate cross-section of typical parameters. But to keep runtimes
// relatively reasonable, I'm mostly not using it.

const NUS: &[f64] = &[3e8, 1e9, 3e9, 3e10, 3e11];
const N_ES: &[f64] = &[1e2, 1e3, 1e4, 1e5, 1e6];
const SS: &[f64] = &[1e0, 1e1, 1e2, 1e3, 1e4];
const THETAS: &[f64] = &[0.05, 0.430, 0.810, 1.190, 1.5707];
const PS: &[f64] = &[1.5, 1.75, 2.5, 3.25, 4.];
const GAMMA_MIN: f64 = 1.;
const GAMMA_MAX: f64 = 1e12;
const GAMMA_CUTOFF: f64 = 1e10;

const LATIN_SQUARE: &[usize] = &[
    1, 4, 2, 3, 0,
    3, 1, 0, 4, 2,
    0, 3, 1, 2, 4,
    2, 0, 4, 1, 3,
    4, 2, 3, 0, 1,
];

fn powerlaw_inner(coeff: Coefficient, stokes: Stokes, row_number: usize) {
    // This function used to loop over every row in the Latin square, but that
    // would lead to the benchmark taking way too long to run.

    let base = row_number * 5;
    let nu = NUS[LATIN_SQUARE[base + 0]];
    let n_e = N_ES[LATIN_SQUARE[base + 1]];
    let s = SS[LATIN_SQUARE[base + 2]];
    let theta = THETAS[LATIN_SQUARE[base + 3]];
    let p = PS[LATIN_SQUARE[base + 4]];

    let bfield = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * nu / (ELECTRON_CHARGE * s);

    rimphony::PowerLawDistribution::new(p)
        .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
        .full_calculation()
        .compute_cgs(coeff, stokes, nu, bfield, n_e, theta);
}


fn powerlaw_ji_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::I, 0);
    });
}

fn powerlaw_jq_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::Q, 0);
    });
}

fn powerlaw_jv_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::V, 0);
    });
}

fn powerlaw_ai_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::I, 0);
    });
}

fn powerlaw_aq_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::Q, 0);
    });
}

fn powerlaw_av_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::V, 0);
    });
}

fn powerlaw_fq_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Faraday, Stokes::Q, 0);
    });
}

fn powerlaw_fv_0(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Faraday, Stokes::V, 0);
    });
}


// A smattering of other parameter configurations.

fn powerlaw_ji_1(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::I, 1);
    });
}

fn powerlaw_jq_2(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::Q, 2);
    });
}

fn powerlaw_jv_3(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Emission, Stokes::V, 3);
    });
}

fn powerlaw_ai_4(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::I, 4);
    });
}

fn powerlaw_aq_1(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::Q, 1);
    });
}

fn powerlaw_av_2(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Absorption, Stokes::V, 2);
    });
}

fn powerlaw_fq_3(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Faraday, Stokes::Q, 3);
    });
}

fn powerlaw_fv_4(b: &mut Bencher) {
    b.iter(|| {
        powerlaw_inner(Coefficient::Faraday, Stokes::V, 4);
    });
}


benchmark_group!(powerlaw_0,
                 powerlaw_ji_0, powerlaw_jq_0, powerlaw_jv_0,
                 powerlaw_ai_0, powerlaw_aq_0, powerlaw_av_0,
                 powerlaw_fq_0, powerlaw_fv_0);
benchmark_group!(powerlaw_others,
                 powerlaw_ji_1, powerlaw_jq_2, powerlaw_jv_3,
                 powerlaw_ai_4, powerlaw_aq_1, powerlaw_av_2,
                 powerlaw_fq_3, powerlaw_fv_4);
benchmark_main!(powerlaw_0, powerlaw_others);
