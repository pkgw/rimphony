// Copyright 2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/// Compute all of the coefficients for the "pitchy kappa" distribution,
/// scaling to cgs values.
extern crate clap;
extern crate rimphony;
extern crate rimphony_test_support;

use rimphony::{Coefficient, Stokes, SynchrotronCalculator};

fn main() {
    let matches = clap::Command::new("all-pitchykappa-cgs")
        .version("0.1.0")
        .about("Compute a set of coefficients for the \"pitchy kappa\" distribution.")
        .arg(
            clap::Arg::new("NU")
                .help("The radio frequency of interest (Hz)")
                .required(true)
                .index(1),
        )
        .arg(
            clap::Arg::new("B")
                .help("The magnetic field (Gauss)")
                .required(true)
                .index(2),
        )
        .arg(
            clap::Arg::new("N_E")
                .help("The electron density (cm^-3)")
                .required(true)
                .index(3),
        )
        .arg(
            clap::Arg::new("THETA")
                .help("The viewing angle relative to the magnetic field (radians)")
                .required(true)
                .index(4),
        )
        .arg(
            clap::Arg::new("KAPPA")
                .help("The kappa parameter of the electron energy distribution")
                .required(true)
                .index(5),
        )
        .arg(
            clap::Arg::new("WIDTH")
                .help("The width parameter of the electron energy distribution")
                .required(true)
                .index(6),
        )
        .arg(
            clap::Arg::new("K")
                .help("The power-law index of the electron pitch-angle distribution")
                .required(true)
                .index(7),
        )
        .get_matches();

    let nu = matches
        .get_one::<String>("NU")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let b = matches
        .get_one::<String>("B")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let n_e = matches
        .get_one::<String>("N_E")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let theta = matches
        .get_one::<String>("THETA")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let kappa = matches
        .get_one::<String>("KAPPA")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let width = matches
        .get_one::<String>("WIDTH")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let k = matches
        .get_one::<String>("K")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    const GAMMA_CUTOFF: f64 = 100.;

    let calc = rimphony::PitchyKappaDistribution::new(kappa, width, k)
        .gamma_cutoff(GAMMA_CUTOFF)
        .full_calculation(rimphony_test_support::default_log());

    println!(
        "    j_I: {:.18e}",
        calc.compute_cgs(Coefficient::Emission, Stokes::I, nu, b, n_e, theta)
    );
    println!(
        "alpha_I: {:.18e}",
        calc.compute_cgs(Coefficient::Absorption, Stokes::I, nu, b, n_e, theta)
    );
    println!(
        "    j_Q: {:.18e}",
        calc.compute_cgs(Coefficient::Emission, Stokes::Q, nu, b, n_e, theta)
    );
    println!(
        "alpha_Q: {:.18e}",
        calc.compute_cgs(Coefficient::Absorption, Stokes::Q, nu, b, n_e, theta)
    );
    println!(
        "    j_V: {:.18e}",
        calc.compute_cgs(Coefficient::Emission, Stokes::V, nu, b, n_e, theta)
    );
    println!(
        "alpha_V: {:.18e}",
        calc.compute_cgs(Coefficient::Absorption, Stokes::V, nu, b, n_e, theta)
    );
    println!(
        "  rho_Q: {:.18e}",
        calc.compute_cgs(Coefficient::Faraday, Stokes::Q, nu, b, n_e, theta)
    );
    println!(
        "  rho_V: {:.18e}",
        calc.compute_cgs(Coefficient::Faraday, Stokes::V, nu, b, n_e, theta)
    );
}
