// Copyright 2017-2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Crank out coefficients for the "pitchy" power law distribution.
extern crate clap;
extern crate rimphony;
extern crate rimphony_test_support;

use rimphony::{PitchyPowerLawDistribution, SynchrotronCalculator};
use rimphony_test_support::Sampler;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let matches = clap::Command::new("crank-out-pitchypl")
        .version("0.1.0")
        .about("Crank out coefficients for random pitchy power law parameters")
        .arg(
            clap::Arg::new("OUTFILE")
                .help("The path of the output file to create")
                .required(true)
                .index(1),
        )
        .arg(
            clap::Arg::new("S_MIN")
                .help("The minimum 's' value to generate")
                .required(true)
                .index(2),
        )
        .arg(
            clap::Arg::new("S_MAX")
                .help("The maximum 's' value to generate")
                .required(true)
                .index(3),
        )
        .arg(
            clap::Arg::new("THETA_MIN")
                .help("The minimum 'theta' value to generate")
                .required(true)
                .index(4),
        )
        .arg(
            clap::Arg::new("THETA_MAX")
                .help("The maximum 'theta' value to generate")
                .required(true)
                .index(5),
        )
        .arg(
            clap::Arg::new("P_MIN")
                .help("The minimum 'p' value to generate")
                .required(true)
                .index(6),
        )
        .arg(
            clap::Arg::new("P_MAX")
                .help("The maximum 'p' value to generate")
                .required(true)
                .index(7),
        )
        .arg(
            clap::Arg::new("K_MIN")
                .help("The minimum 'k' value to generate")
                .required(true)
                .index(8),
        )
        .arg(
            clap::Arg::new("K_MAX")
                .help("The maximum 'k' value to generate")
                .required(true)
                .index(9),
        )
        .get_matches();

    let outfile = matches.get_one::<PathBuf>("OUTFILE").unwrap();
    let log = rimphony_test_support::default_log();

    let s_sampler = Sampler::new(
        true,
        matches
            .get_one::<String>("S_MIN")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
        matches
            .get_one::<String>("S_MAX")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
    );
    let theta_sampler = Sampler::new(
        false,
        matches
            .get_one::<String>("THETA_MIN")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
        matches
            .get_one::<String>("THETA_MAX")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
    );
    let p_sampler = Sampler::new(
        false,
        matches
            .get_one::<String>("P_MIN")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
        matches
            .get_one::<String>("P_MAX")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
    );
    let k_sampler = Sampler::new(
        false,
        matches
            .get_one::<String>("K_MIN")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
        matches
            .get_one::<String>("K_MAX")
            .unwrap()
            .parse::<f64>()
            .unwrap(),
    );

    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .append(true)
        .open(outfile)
        .unwrap();

    writeln!(
        file,
        "s(log)\t\
                    theta(lin)\t\
                    p(lin)\t\
                    k(lin)\t\
                    time_ms(meta)\t\
                    j_I(res)\t\
                    alpha_I(res)\t\
                    j_Q(res)\t\
                    alpha_Q(res)\t\
                    j_V(res)\t\
                    alpha_V(res)\t\
                    rho_Q(res)\t\
                    rho_V(res)"
    )
    .expect("write error");

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
            .full_calculation(log.clone())
            .compute_all_dimensionless(s, theta);
        let elapsed = t0.elapsed();
        let ms = elapsed.as_secs() as f64 * 1000. + elapsed.subsec_nanos() as f64 * 1e-6;

        writeln!(
            file,
            "{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t\
                  {:.16e}",
            s,
            theta,
            p,
            k,
            ms,
            vals[0],
            vals[1],
            vals[2],
            vals[3],
            vals[4],
            vals[5],
            vals[6],
            vals[7]
        )
        .expect("write error");
    }
}
