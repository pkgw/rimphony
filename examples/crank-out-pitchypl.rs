// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Crank out coefficients for the "pitchy" power law distribution.

extern crate clap;
extern crate rand;
extern crate rimphony;

use rimphony::{PitchyPowerLawDistribution, SynchrotronCalculator};
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

struct Sampler {
    is_log: bool,
    low: f64,
    range: f64
}

impl Sampler {
    fn new(is_log: bool, mut low: f64, mut high: f64) -> Self {
        if low > high {
            let tmp = high;
            high = low;
            low = tmp;
        }

        if is_log {
            low = low.ln();
            high = high.ln();
        }

        Sampler { is_log: is_log, low: low, range: high - low }
    }

    fn get(&self) -> f64 {
        let n = self.low + rand::random::<f64>() * self.range;

        if self.is_log {
            n.exp()
        } else {
            n
        }
    }
}


fn main() {
    let matches = clap::App::new("crank-out-pitchypl")
        .version("0.1.0")
        .about("Crank out coefficients for random pitchy power law parameters")
        .arg(clap::Arg::with_name("OUTFILE")
             .help("The path of the output file to create")
             .required(true)
             .index(1))
        .arg(clap::Arg::with_name("S_MIN")
             .help("The minimum 's' value to generate")
             .required(true)
             .index(2))
        .arg(clap::Arg::with_name("S_MAX")
             .help("The maximum 's' value to generate")
             .required(true)
             .index(3))
        .arg(clap::Arg::with_name("THETA_MIN")
             .help("The minimum 'theta' value to generate")
             .required(true)
             .index(4))
        .arg(clap::Arg::with_name("THETA_MAX")
             .help("The maximum 'theta' value to generate")
             .required(true)
             .index(5))
        .arg(clap::Arg::with_name("P_MIN")
             .help("The minimum 'p' value to generate")
             .required(true)
             .index(6))
        .arg(clap::Arg::with_name("P_MAX")
             .help("The maximum 'p' value to generate")
             .required(true)
             .index(7))
        .arg(clap::Arg::with_name("K_MIN")
             .help("The minimum 'k' value to generate")
             .required(true)
             .index(8))
        .arg(clap::Arg::with_name("K_MAX")
             .help("The maximum 'k' value to generate")
             .required(true)
             .index(9))
        .get_matches();

    let outfile = PathBuf::from(matches.value_of_os("OUTFILE").unwrap());

    let s_sampler = Sampler::new(
        true,
        matches.value_of("S_MIN").unwrap().parse::<f64>().unwrap(),
        matches.value_of("S_MAX").unwrap().parse::<f64>().unwrap());
    let theta_sampler = Sampler::new(
        false,
        matches.value_of("THETA_MIN").unwrap().parse::<f64>().unwrap(),
        matches.value_of("THETA_MAX").unwrap().parse::<f64>().unwrap());
    let p_sampler = Sampler::new(
        false,
        matches.value_of("P_MIN").unwrap().parse::<f64>().unwrap(),
        matches.value_of("P_MAX").unwrap().parse::<f64>().unwrap());
    let k_sampler = Sampler::new(
        false,
        matches.value_of("K_MIN").unwrap().parse::<f64>().unwrap(),
        matches.value_of("K_MAX").unwrap().parse::<f64>().unwrap());

    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .append(true)
        .open(outfile)
        .unwrap();

    writeln!(file, "# s(log) theta(lin) p(lin) k(lin) !time_ms").expect("write error");

    loop {
        let s = s_sampler.get();
        let theta = theta_sampler.get();
        let p = p_sampler.get();
        let k = k_sampler.get();

        const GAMMA_MIN: f64 = 1.;
        const GAMMA_MAX: f64 = 1e12;
        const GAMMA_CUTOFF: f64 = 1e10;

        let t0 = Instant::now();
        let vals = PitchyPowerLawDistribution::new(p, k)
            .gamma_limits(GAMMA_MIN, GAMMA_MAX, GAMMA_CUTOFF)
            .full_calculation()
            .compute_all_dimensionless(s, theta);
        let elapsed = t0.elapsed();
        let ms = elapsed.as_secs() as f64 * 1000. + elapsed.subsec_nanos() as f64 * 1e-6;

        writeln!(file,
                 "{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}",
                 s, theta, p, k, ms,
                 vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]
        ).expect("write error");
    }
}
