// Copyright 2017-2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Check our numbers against Symphony's.
///
/// This test suite compares against a random subset of a larger database of
/// examples computed with Symphony; this is a compromise between coverage
/// and runtime.

extern crate rand;
extern crate regex;
extern crate rimphony;
extern crate rimphony_test_support;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI,
               Coefficient, Stokes, SynchrotronCalculator};
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

const TOP: &'static str = env!("CARGO_MANIFEST_DIR");

/// These computations were parameterized by the harmonic number *s*, observer
/// angle *theta*, and power-law index *p*. Other parameters were chosen
/// fiducially; the coefficients can be scaled deterministically to apply for
/// any values of *B*, *nu*, and *n_e*.

fn compare_powerlaw_subset(coeff: Coefficient, stokes: Stokes, cname: &str, index: usize, fraction_to_check: f64, rtol: f64) {
    const NU: f64 = 1e9;
    const N_E: f64 = 1.;
    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1e12;
    const GAMMA_CUTOFF: f64 = 1e10;

    let log = rimphony_test_support::default_log();
    let re = regex::Regex::new(r"\s+").unwrap();

    let mut p = PathBuf::from(TOP);
    p.push("tests");
    p.push("symphony-powerlaw.txt");

    let f = BufReader::new(File::open(p).unwrap());

    for line in f.lines() {
        if rand::random::<f64>() > fraction_to_check {
            continue;
        }

        let v: Vec<f64> = re.split(line.unwrap().as_ref()).map(|s| s.parse::<f64>().unwrap()).collect();
        let s = v[0];
        let theta = v[1];
        let p = v[2];
        let bfield = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * s);

        let t0 = Instant::now();
        let ours = rimphony::PowerLawDistribution::new(p)
            .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
            .full_calculation(log.clone())
            .compute_cgs(coeff, stokes, NU, bfield, N_E, theta);

        let elapsed = t0.elapsed();
        if elapsed.as_secs() > 3 {
            println!("SLOW: c = {:?}/{:?} s = {:.10e} theta = {:.10e} p = {:.10e} elapsed = {:.0} ms",
                     coeff, stokes, s, theta, p,
                     elapsed.as_secs() as f64 * 1000. + elapsed.subsec_nanos() as f64 * 1e-6);
        }

        let theirs = v[index];
        let rel_err = ((ours - theirs) / theirs).abs();

        if rel_err > rtol {
            panic!("disagree with Symphony for {}; they have {:.6e} , we have {:.6e}; at parameters:

s: {:.10e}
theta: {:.10e}
p: {:.10e}", cname, theirs, ours, s, theta, p);
        }
    }
}

const TOL: f64 = 0.01;

#[test]
fn compare_powerlaw_subset_ji() {
    compare_powerlaw_subset(Coefficient::Emission, Stokes::I, "J_I", 3, 0.03, TOL);
}

#[test]
fn compare_powerlaw_subset_jq() {
    compare_powerlaw_subset(Coefficient::Emission, Stokes::Q, "J_Q", 5, 0.03, TOL);
}

#[test]
fn compare_powerlaw_subset_jv() {
    compare_powerlaw_subset(Coefficient::Emission, Stokes::V, "J_V", 7, 0.03, TOL);
}

#[test]
fn compare_powerlaw_subset_ai() {
    compare_powerlaw_subset(Coefficient::Absorption, Stokes::I, "A_I", 4, 0.03, TOL);
}

#[test]
fn compare_powerlaw_subset_aq() {
    compare_powerlaw_subset(Coefficient::Absorption, Stokes::Q, "A_Q", 6, 0.03, TOL);
}

#[test]
fn compare_powerlaw_subset_av() {
    compare_powerlaw_subset(Coefficient::Absorption, Stokes::V, "A_V", 8, 0.03, TOL);
}
