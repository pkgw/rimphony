// Copyright 2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Compute isotropic power-law coefficients for some canned, smoothly-varying
/// input parameters.
///
/// The output of this tool can be fed into neurosynchro’s testing framework
/// to see how well its approximations perform in a more or less end-to-end
/// test.

#[macro_use]
extern crate clap;
extern crate rimphony;
extern crate rimphony_test_support;
extern crate slog;

use rimphony::{PowerLawDistribution, SynchrotronCalculator};
use slog::Logger;
use std::time::Instant;

const GAMMA_MIN: f64 = 1.;
const GAMMA_MAX: f64 = 1e12;
const GAMMA_CUTOFF: f64 = 1e10;

#[derive(Debug)]
struct Record {
    d: f64,
    s: f64,
    theta: f64,
    psi: f64,
    n_e: f64,

    time_ms: f64,
    coeffs: [f64; 8],

    // Specific to this distribution:
    p: f64,
}

impl Record {
    pub fn new(d: f64, s: f64, theta: f64, psi: f64, n_e: f64, p: f64, log: &mut Logger) -> Self {
        let t0 = Instant::now();
        let coeffs = PowerLawDistribution::new(p)
            .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
            .full_calculation(log.clone())
            .compute_all_dimensionless(s, theta);
        let elapsed = t0.elapsed();
        let time_ms = elapsed.as_secs() as f64 * 1000. + elapsed.subsec_nanos() as f64 * 1e-6;

        Record {
            d,
            s,
            theta,
            psi,
            n_e,
            p,
            time_ms,
            coeffs,
        }
    }
}

/// This demo prints out 64 sets of coefficients, with gentle linear
/// interpolation over a pretty small range of parameter values.
#[derive(Debug)]
struct AlmostUniform1Demo {
    n: usize,
    log: Logger,
}

impl AlmostUniform1Demo {
    pub fn new(log: Logger) -> Self {
        AlmostUniform1Demo { n: 0, log }
    }
}

impl Iterator for AlmostUniform1Demo {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        const STEPS: usize = 64;

        if self.n >= STEPS {
            return None;
        }

        let x = self.n as f64 / (STEPS - 1) as f64; // => [0, 1]

        let d = x * 3e10; // centimeters; R_jup ~= 7e9 cm
        let s = 100. - 10. * x;
        let theta = 0.5 + 0.1 * x;
        let psi = 0. + 0.1 * x;
        let n_e = 1e5 - 3e4 * x;
        let p = 3. - 0.5 * x;

        self.n += 1;
        Some(Record::new(d, s, theta, psi, n_e, p, &mut self.log))
    }
}

fn main() {
    let matches = clap::App::new(crate_name!())
        .version(crate_version!())
        .about("Compute coefficients for some canned test cases")
        .arg(
            clap::Arg::with_name("DEMONAME")
                .help("Which demo to compute")
                .required(true)
                .possible_values(&["almostuniform1"])
                .index(1),
        )
        .get_matches();

    let log = rimphony_test_support::default_log();

    let gen: Box<Iterator<Item = Record>> = Box::new(match matches.value_of("DEMONAME").unwrap() {
        "almostuniform1" => AlmostUniform1Demo::new(log),
        _ => unreachable!(),
    });

    println!(
        "s(lin)\t\
              theta(lin)\t\
              p(lin)\t\
              d(meta)\t\
              psi(meta)\t\
              n_e(meta)\t\
              time_ms(meta)\t\
              j_I(res)\t\
              alpha_I(res)\t\
              j_Q(res)\t\
              alpha_Q(res)\t\
              j_V(res)\t\
              alpha_V(res)\t\
              rho_Q(res)\t\
              rho_V(res)"
    );

    for rec in gen {
        let c = &rec.coeffs;
        println!(
            "{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}\t{:.16e}\t{:.16e}",
            rec.s,
            rec.theta,
            rec.p,
            rec.d,
            rec.psi,
            rec.n_e,
            rec.time_ms,
            c[0],
            c[1],
            c[2],
            c[3],
            c[4],
            c[5],
            c[6],
            c[7]
        );
    }
}
