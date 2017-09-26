// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! Calculate radiative transfer coefficients for synchrotron emission.

This crate is a Rust port of
[symphony](https://github.com/AFD-Illinois/symphony), a synchrotron code
written in C that is primarily the work of Alex Pandya and the group of
Charles Gammie at the University of Illinois. This port attempts to clean
up the code to make it a bit easier to experiment with. The key
publications are [Pandya, Zhang, Chandra, and Gammie (2016;
DOI:10.3847/0004-637X/822/1/34](https://dx.doi.org/10.3847/0004-637X/822/1/34)
and [Leung, Gammie, and Noble
(2011)](https://dx.doi.org/10.1088/0004-637X/737/1/21).

The basic structure of the problem is that we need to do an integral in a 2D
quarter-plane defined by the variables *gamma* (>= 1) and *n* (>= 1).
Technically, *n* can only take on discrete values, but we generally have to
integrate over very large values of *n* such that we can treat it as being
continous.

*/

#![deny(missing_docs)]

extern crate gsl_sys;
extern crate leung_bessel;

#[cfg(test)] #[macro_use] extern crate assert_approx_eq;
#[cfg(test)] extern crate rand;

use std::f64;

mod gsl;
mod symphony;

pub use f64::consts::PI;

/// Two times pi, as an `f64`.
pub const TWO_PI: f64 = 2. * PI;

/// The mass of the electron in cgs (grams).
pub const MASS_ELECTRON: f64 = 9.1093826e-28;

/// The speed of light in cgs (centimeters per second).
pub const SPEED_LIGHT: f64 = 2.99792458e10;

/// The charge of the electron, in cgs (esu's).
pub const ELECTRON_CHARGE: f64 = 4.80320680e-10;


/// Which Stokes parameter we are analyzing. There's no "U" option because our
/// linear polarization basis is defined such that all U terms are zero; one
/// can fairly easily rotate these parameters into a fixed frame if the
/// magnetic field rotates over the line of sight.
#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub enum Stokes {
    /// The total intensity.
    I,

    /// The linear polarization.
    Q,

    /// The circular polarization. The sign convention aspires to agree with
    /// the IEEE/IAU convention, where `V = RCP - LCP` and RCP means that the
    /// electric field rotates clockwise as seen from the transmitter. Sign
    /// errors in Stokes V are endemic.
    V
}


/// Which radiative transfer coefficient to compute.
#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub enum Coefficient {
    /// The emission cofficient, in units of ergs per second per square
    /// centimeter per Hertz per steradian.
    Emission,

    /// The absorption coefficient, in units of inverse centimeters.
    Absorption,
}


/// An electron distribution function. We can compute synchrotron coefficients
/// using different distributions.
pub trait DistributionFunction {
    /// This function gives the modified distribution function `d(n_e) /
    /// (gamma^2 beta d(gamma) d(cos xi) d(phi)])`. In all the cases in
    /// Symphony, isotropy and gyrotropy are assumed, so the `d(cos xi)
    /// d(phi)` become 1/4pi. The funny scalings are due to me being stuck on
    /// the numerical derivative step. To be revisited.
    fn calc_f(&self, gamma: f64, cos_xi: f64) -> f64;

    /// The derivative of `calc_f` with regards to `gamma` and `cos_xi`. This
    /// must include the derivative with regards to the `1 / (gamma^2 beta)`
    /// factor that turns the gamma/cos xi/phi coordinates into p^3
    /// coordinates.
    fn calc_f_derivatives(&self, gamma: f64, cos_xi: f64) -> (f64, f64);
}


/// A type that can compute synchrotron coefficients.
pub trait SynchrotronCalculator {
    /// Compute a coefficient using the dimensionless parameters that
    /// characterize the problem. *s* is the harmonic number of the observing
    /// frequency relative to the cyclotron frequency and *theta* is the
    /// observer angle relative to the magnetic field direction (in radians).
    fn compute_dimensionless(&self, coeff: Coefficient, stokes: Stokes,
                             s: f64, theta: f64) -> f64;

    /// Compute a coefficient using standard physics parameters in cgs units.
    /// *nu* is the observing frequency in Hz, *b* is the magnetic field
    /// strength in Gauss, *n_e* is the electron density in cm^-3, and *theta*
    /// is the observer angle relative to the magnetic field direction in
    /// radians.
    fn compute_cgs(&self, coeff: Coefficient, stokes: Stokes,
                   nu: f64, b: f64, n_e: f64, theta: f64) -> f64 {
        let nu_c = ELECTRON_CHARGE * b / (TWO_PI * MASS_ELECTRON * SPEED_LIGHT);
        let val = self.compute_dimensionless(coeff, stokes, nu / nu_c, theta);

        match coeff {
            Coefficient::Emission => val * n_e * nu,
            Coefficient::Absorption => val * n_e / nu,
        }
    }

    /// Compute all six coefficients in their dimensionless form. The return
    /// value is an array of `[j_i, alpha_i, j_q, alpha_q, j_v, alpha_v]`.
    fn compute_all_dimensionless(&self, s: f64, theta: f64) -> [f64; 6] {
        let mut rv = [0_f64; 6];

        rv[0] = self.compute_dimensionless(Coefficient::Emission, Stokes::I, s, theta);
        rv[1] = self.compute_dimensionless(Coefficient::Absorption, Stokes::I, s, theta);
        rv[2] = self.compute_dimensionless(Coefficient::Emission, Stokes::Q, s, theta);
        rv[3] = self.compute_dimensionless(Coefficient::Absorption, Stokes::Q, s, theta);
        rv[4] = self.compute_dimensionless(Coefficient::Emission, Stokes::V, s, theta);
        rv[5] = self.compute_dimensionless(Coefficient::Absorption, Stokes::V, s, theta);

        rv
    }

    /// Compute all six coefficients in their cgs form. The return
    /// value is an array of `[j_i, alpha_i, j_q, alpha_q, j_v, alpha_v]`.
    fn compute_all_cgs(&self, nu: f64, b: f64, n_e: f64, theta: f64) -> [f64; 6] {
        let mut rv = [0_f64; 6];

        rv[0] = self.compute_cgs(Coefficient::Emission, Stokes::I, nu, b, n_e,  theta);
        rv[1] = self.compute_cgs(Coefficient::Absorption, Stokes::I, nu, b, n_e,  theta);
        rv[2] = self.compute_cgs(Coefficient::Emission, Stokes::Q, nu, b, n_e,  theta);
        rv[3] = self.compute_cgs(Coefficient::Absorption, Stokes::Q, nu, b, n_e,  theta);
        rv[4] = self.compute_cgs(Coefficient::Emission, Stokes::V, nu, b, n_e,  theta);
        rv[5] = self.compute_cgs(Coefficient::Absorption, Stokes::V, nu, b, n_e,  theta);

        rv
    }
}


// Distributions

pub mod power_law;
pub use power_law::PowerLawDistribution;

pub mod pitchy_pl;
pub use pitchy_pl::PitchyPowerLawDistribution;


/// The FullSunchrotronCalculator implements the fully detailed
/// double-integral calculation.
#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub struct FullSynchrotronCalculator<D>(D);


impl<D: DistributionFunction> SynchrotronCalculator for FullSynchrotronCalculator<D> {
    fn compute_dimensionless(&self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> f64 {
        symphony::SymphonyCalculationState::new(&self.0, coeff, stokes, s, theta)
            .compute()
    }
}
