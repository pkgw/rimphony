/// Check our numbers against Symphony's.

extern crate rand;
extern crate rimphony;

use rimphony::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI, Coefficient, Stokes};
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::path::PathBuf;
//use std::str::StrExt;

const TOP: &'static str = env!("CARGO_MANIFEST_DIR");

/// These computations were parameterized by the harmonic number *s*, observer
/// angle *theta*, and power-law index *p*. Other parameters were chosen
/// fiducially; the coefficients can be scaled deterministically to apply for
/// any values of *B*, *nu*, and *n_e*.

fn compare_powerlaw_subset(coeff: rimphony::Coefficient, cname: &str, index: usize, fraction_to_check: f64, rtol: f64) {
    const NU: f64 = 1e9;
    const N_E: f64 = 1.;
    const GAMMA_MIN: f64 = 1.;
    const GAMMA_MAX: f64 = 1000.;
    const GAMMA_CUTOFF: f64 = 1e7;

    let mut p = PathBuf::from(TOP);
    p.push("tests");
    p.push("symphony-powerlaw.txt");

    let f = BufReader::new(File::open(p).unwrap());

    for line in f.lines() {
        if rand::random::<f64>() > fraction_to_check {
            continue;
        }

        let v: Vec<f64> = line.unwrap().split("\t").map(|s| s.parse::<f64>().unwrap()).collect();
        let s = v[0];
        let theta = v[1];
        let p = v[2];
        let bfield = TWO_PI * MASS_ELECTRON * SPEED_LIGHT * NU / (ELECTRON_CHARGE * s);

        let ours = rimphony::PowerLawDistribution::new(p)
            .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
            .finish(coeff, NU, bfield, N_E, theta)
            .compute();
        let theirs = v[index];
        let rel_err = ((ours - theirs) / theirs).abs();

        if rel_err > rtol {
            panic!("disagree with Symphony for {}; they have {:.6e}, we have {:.6e}; at parameters:

s: {:.6e}
theta: {:.6e}
p: {:.6e}", cname, theirs, ours, s, theta, p);
        }
    }
}

#[test]
fn compare_powerlaw_subset_ji() {
    compare_powerlaw_subset(Coefficient::Emission(Stokes::I), "J_I", 3, 0.002, 1e-4);
}

#[test]
fn compare_powerlaw_subset_jq() {
    compare_powerlaw_subset(Coefficient::Emission(Stokes::Q), "J_Q", 5, 0.002, 1e-4);
}

/// Stokes V is a bit more demanding so the error threshold is looser.
#[test]
fn compare_powerlaw_subset_jv() {
    compare_powerlaw_subset(Coefficient::Emission(Stokes::V), "J_V", 7, 0.002, 1e-3);
}

#[test]
fn compare_powerlaw_subset_ai() {
    compare_powerlaw_subset(Coefficient::Absorption(Stokes::I), "A_I", 4, 0.002, 1e-4);
}

#[test]
fn compare_powerlaw_subset_aq() {
    compare_powerlaw_subset(Coefficient::Absorption(Stokes::Q), "A_Q", 6, 0.002, 1e-4);
}

/// Stokes V is a bit more demanding so the error threshold is looser.
#[test]
fn compare_powerlaw_subset_av() {
    compare_powerlaw_subset(Coefficient::Absorption(Stokes::V), "A_V", 8, 0.002, 1e-3);
}
