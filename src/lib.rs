// Copyright 2017-2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! Calculate radiative transfer coefficients for synchrotron emission.

This crate is a Rust port of
[symphony](https://github.com/AFD-Illinois/symphony), a synchrotron code
written in C that is primarily the work of Alex Pandya and the group of
Charles Gammie at the University of Illinois. This port attempts to clean
up the code to make it a bit easier to experiment with. The key
publications are [Pandya, Zhang, Chandra, and Gammie (2016;
DOI:10.3847/0004-637X/822/1/34)](https://dx.doi.org/10.3847/0004-637X/822/1/34)
and [Leung, Gammie, and Noble
(2011; DOI:10.1088/0004-637X/737/1/21)](https://dx.doi.org/10.1088/0004-637X/737/1/21).

It has been expanded to compute Faraday rotation and conversion coefficients
using the formalism developed by [Heyvaerts et al. (2013;
DOI:10.1093/mnras/stt135)](https://dx.doi.org/10.1093/mnras/stt135).

The basic structure of the problem is that we need to do a 2D integral. In
Symphony, this integral is over a quarter-plane defined by the variables γ (>=
1) and *n* (>= 1). Technically, *n* can only take on discrete values, but we
generally have to integrate over very large values of *n* such that we can
treat it as being continous. In the Heyvaerts formalism, the integral is over
a parabolic surface in two variables called σ and ϖ.

To calculate coefficients, create an instance of the one of the “distribution”
structs defined below, then use one of the methods on the struct to obtain a
helper object implementing the
[SynchrotronCalculator](trait.SynchrotronCalculator.html) trait. In some
cases, you can obtain different kinds of calculators: ones that do the fully
detailed calculations, or ones that use various approximations that can help
with checking results.

*/

#![deny(missing_docs)]

extern crate gsl_sys;
extern crate leung_bessel;
extern crate special_fun;
#[macro_use] extern crate slog;

#[cfg(test)] #[macro_use] extern crate assert_approx_eq;
#[cfg(test)] extern crate rand;
#[cfg(test)] extern crate rimphony_test_support;

use std::f64;
use std::fmt;

mod gsl;
mod heyvaerts;
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

    /// A Faraday conversion or rotation coefficient, in units of inverse
    /// centimeters. The "Stokes Q" coefficient is the Faraday conversion
    /// coefficient, sometimes denoted `rho_{UV}`, and the "Stokes V"
    /// coefficient is the Faraday rotation coefficient, sometimes denoted
    /// `rho_{QU}`. It is not meaningful to refer to the "Stokes I" Faraday
    /// coefficient.
    Faraday,
}

/// An electron distribution function. We can compute synchrotron coefficients
/// using different distributions.
pub trait DistributionFunction: fmt::Debug {
    /// This function gives the modified distribution function `d(n_e) /
    /// [d(gamma) d(cos xi) d(phi)]`. In all the cases in Symphony, isotropy
    /// and gyrotropy are assumed, so the `d(cos xi) d(phi)` becomes 1/4π.
    /// This function should be normalized such that its integral over all of
    /// momentum space is 1/(mc)³. In the gamma/cos xi/phi parameterization,
    /// that integral is:
    ///
    /// ```text
    /// (mc)^{-3} \int_0^{2\pi} d\phi \int_{-1}^{1} d\cos\xi \int_1^\infty d\gamma \gamma^2 \beta f
    /// ```
    ///
    /// So in the isotropic case the requirement is that
    ///
    /// ```text
    /// \int_1^\infty d\gamma \gamma^2 \beta f = 1 / 4 \pi
    /// ```
    ///
    /// (The factor of (mc)³ is simplifying the fact that the distribution
    /// function should really have units of inverse momentum cubed, since the
    /// best definition of the distribution function is as the probability
    /// density function that a given particle has a particular location in
    /// momentum space:
    ///
    /// ```text
    /// \int d^3p f = 1
    /// ```
    ///
    /// In practice virtually everyone expresses their distribution functions
    /// in terms of γ, so the (mc)³ factor is pretty much universal.)
    fn calc_f(&self, gamma: f64, cos_xi: f64) -> f64;

    /// The partial derivatives of `calc_f` with regards to `gamma` and
    /// `cos_xi`.
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
            Coefficient::Faraday => val * n_e / nu,
        }
    }

    /// Compute all eight coefficients in their dimensionless form. The return
    /// value is an array of `[j_i, alpha_i, j_q, alpha_q, j_v, alpha_v,
    /// faraday_q, faraday_v]`.
    fn compute_all_dimensionless(&self, s: f64, theta: f64) -> [f64; 8] {
        let mut rv = [0_f64; 8];

        rv[0] = self.compute_dimensionless(Coefficient::Emission, Stokes::I, s, theta);
        rv[1] = self.compute_dimensionless(Coefficient::Absorption, Stokes::I, s, theta);
        rv[2] = self.compute_dimensionless(Coefficient::Emission, Stokes::Q, s, theta);
        rv[3] = self.compute_dimensionless(Coefficient::Absorption, Stokes::Q, s, theta);
        rv[4] = self.compute_dimensionless(Coefficient::Emission, Stokes::V, s, theta);
        rv[5] = self.compute_dimensionless(Coefficient::Absorption, Stokes::V, s, theta);
        rv[6] = self.compute_dimensionless(Coefficient::Faraday, Stokes::Q, s, theta);
        rv[7] = self.compute_dimensionless(Coefficient::Faraday, Stokes::V, s, theta);

        rv
    }

    /// Compute all six coefficients in their cgs form. The return
    /// value is an array of `[j_i, alpha_i, j_q, alpha_q, j_v, alpha_v,
    /// faraday_q, faraday_v]`.
    fn compute_all_cgs(&self, nu: f64, b: f64, n_e: f64, theta: f64) -> [f64; 8] {
        let mut rv = [0_f64; 8];

        rv[0] = self.compute_cgs(Coefficient::Emission, Stokes::I, nu, b, n_e, theta);
        rv[1] = self.compute_cgs(Coefficient::Absorption, Stokes::I, nu, b, n_e, theta);
        rv[2] = self.compute_cgs(Coefficient::Emission, Stokes::Q, nu, b, n_e, theta);
        rv[3] = self.compute_cgs(Coefficient::Absorption, Stokes::Q, nu, b, n_e, theta);
        rv[4] = self.compute_cgs(Coefficient::Emission, Stokes::V, nu, b, n_e, theta);
        rv[5] = self.compute_cgs(Coefficient::Absorption, Stokes::V, nu, b, n_e, theta);
        rv[6] = self.compute_cgs(Coefficient::Faraday, Stokes::Q, nu, b, n_e, theta);
        rv[7] = self.compute_cgs(Coefficient::Faraday, Stokes::V, nu, b, n_e, theta);

        rv
    }
}


// Distributions

pub mod power_law;
pub use power_law::PowerLawDistribution;

pub mod pitchy_pl;
pub use pitchy_pl::PitchyPowerLawDistribution;

pub mod pitchy_kappa;
pub use pitchy_kappa::PitchyKappaDistribution;

pub mod thermal_juettner;
pub use thermal_juettner::ThermalJuettnerDistribution;


/// The FullSunchrotronCalculator implements the fully detailed
/// double-integral calculation.
#[derive(Clone, Debug)]
pub struct FullSynchrotronCalculator<D> {
    distrib: D,
    logger: slog::Logger,
}

impl<D: DistributionFunction> SynchrotronCalculator for FullSynchrotronCalculator<D> {
    fn compute_dimensionless(&self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> f64 {
        match (coeff, stokes) {
            (Coefficient::Faraday, Stokes::I) =>
                f64::NAN,
            (Coefficient::Emission, _)|(Coefficient::Absorption, _) =>
                symphony::compute_dimensionless(&self.distrib, &self.logger, coeff, stokes, s, theta),
            (Coefficient::Faraday, _) =>
                heyvaerts::compute_dimensionless(&self.distrib, coeff, stokes, s, theta), // XXX LOG
        }
    }
}

impl<D: DistributionFunction> FullSynchrotronCalculator<D> {
    /// A diagnostic function for investigating the Symphony N integral.
    ///
    /// This function is only of interest if you're examining the internal
    /// behavior of the code.
    pub fn diagnostic_symphony_n_integral(
        &self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64, n_lo: f64, n_hi: f64
    ) -> gsl::GslResult<f64> {
        symphony::diagnostic_n_integral(
            &self.distrib, &self.logger, coeff, stokes, s, theta, n_lo, n_hi
        )
    }

    /// A diagnostic function for investigating the Symphony gamma integral.
    ///
    /// This function is only of interest if you're examining the internal
    /// behavior of the code.
    pub fn diagnostic_symphony_gamma_integral(
        &self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64, n: f64
    ) -> f64 {
        symphony::diagnostic_gamma_integral(
            &self.distrib, &self.logger, coeff, stokes, s, theta, n
        )
    }

    /// A diagnostic function for investigating the Symphony gamma integral.
    ///
    /// This function is only of interest if you're examining the internal
    /// behavior of the code.
    pub fn diagnostic_symphony_gamma_integrand(
        &self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64, n: f64, gamma: f64
    ) -> f64 {
        symphony::diagnostic_gamma_integrand(
            &self.distrib, &self.logger, coeff, stokes, s, theta, n, gamma
        )
    }

    /// A diagnostic function for investigating the Symphony double integral.
    ///
    /// This reverses the order of the Symphony integral, and returns the
    /// contribution across all relevant *n* at fixed gamma. In the Symphony
    /// approach, the final coefficient is integrated across all relevant
    /// gamma at fixed *n*.
    pub fn diagnostic_symphony_gamma_contribution(
        &self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64, gamma: f64
    ) -> f64 {
        symphony::diagnostic_gamma_contribution(
            &self.distrib, &self.logger, coeff, stokes, s, theta, gamma
        )
    }
}
