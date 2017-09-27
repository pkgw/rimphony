// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! A thermal Jüttner distribution.

The electrons are isotropic. *T* is the ratio of the particles’ thermal energy
to their rest-mass energy.

*/

use special_fun::FloatSpecial;
use std::f64;

use gsl;
use super::{Coefficient, DistributionFunction, FullSynchrotronCalculator, Stokes, SynchrotronCalculator};
use super::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI};


/// Parameters for a thermal Jüttner electron distribution. See the module-level
/// documentation for details.
#[derive(Copy,Clone,Debug,PartialEq)]
pub struct ThermalJuettnerDistribution {
    neg_inverse_t: f64,
    norm: f64,
}


impl DistributionFunction for ThermalJuettnerDistribution {
    fn calc_f(&self, gamma: f64, _cos_xi: f64) -> f64 {
        self.norm * (self.neg_inverse_t * gamma).exp()
    }

    fn calc_f_derivatives(&self, gamma: f64, _cos_xi: f64) -> (f64, f64) {
        let dfdg = self.norm * (self.neg_inverse_t * gamma).exp() * self.neg_inverse_t;
        let dfdcx = 0.;
        (dfdg, dfdcx)
    }
}


impl ThermalJuettnerDistribution {
    /// Create a new thermal Jüttner distribution with the specified
    /// dimensionless temperature.
    pub fn new(t: f64) -> Self {
        ThermalJuettnerDistribution {
            neg_inverse_t: -1. / t,
            norm: f64::NAN,
        }
    }

    /// Compute the normalization factor. The definition of our distribution
    /// function is such that the integral over gamma, with the momentum-space
    /// conversion factor of gamma^2*beta, should be 1 / 4pi (since we’re
    /// isotropic).
    fn normalize(&mut self) {
        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let integral = ws.qagiu(|g| g * (g.powi(2) - 1.).sqrt() * (self.neg_inverse_t * g).exp(), 1.)
            .tolerance(0., 1e-5)
            .compute()
            .map(|r| r.value)
            .unwrap();
        self.norm = 1. / (2. * TWO_PI * integral);
    }

    /// Create a SynchrotronCalculator from this set of parameters. The
    /// calculator will use the full, detailed double integral calculation to
    /// evaluate all coefficients.
    pub fn full_calculation(mut self) -> FullSynchrotronCalculator<Self> {
        self.normalize();
        FullSynchrotronCalculator(self)
    }

    /// Create a SynchrotronCalculator from this set of parameters that uses
    /// high-frequency approximations to compute coefficients quickly. It is
    /// *not* verified that the input parameters are necessarily in the bounds
    /// for which the approximation is accurate.
    pub fn high_freq_approximation(&mut self) -> HighFrequencyApproximation {
        self.normalize();
        HighFrequencyApproximation(self)
    }
}


/// This struct allows one to compute synchrotron coefficients using the
/// approximations that are appropriate for power-law distributions at high
/// resonance numbers.
///
/// Note that no bounds-checking is performed to validate that the input
/// parameters are ones such that the approximations are good.
#[derive(Copy,Clone,Debug,PartialEq)]
pub struct HighFrequencyApproximation<'a>(&'a ThermalJuettnerDistribution);

impl<'a> SynchrotronCalculator for HighFrequencyApproximation<'a> {
    fn compute_dimensionless(&self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> f64 {
        match (coeff, stokes) {
            (Coefficient::Faraday, Stokes::I) => f64::NAN,
            (Coefficient::Faraday, Stokes::Q) => self.faraday_q(s, theta),
            (Coefficient::Faraday, Stokes::V) => self.faraday_v(s, theta),
            (_, _) => f64::NAN,
        }
    }
}

impl<'a> HighFrequencyApproximation<'a> {
    /// From Heyvaerts+ (2013) equation 43.
    ///
    /// For this approximation to work, *s* must be sufficiently large to make
    /// the HF approximation valid, *and* T must be above 10 or so, as per the
    /// left panel of Heyvaerts Figure 2. Larger temperatures lead to much
    /// larger minimum values of *s* such that values of >~ 5e4 are necessary.
    fn faraday_q(&self, s: f64, theta: f64) -> f64 {
        // In the Heyvaerts equation there's a ratio of `omega_plasma^2 /
        // omega`; our calculations are framed in terms of `n_e` and `nu`,
        // which means that the prefactor becomes:
        const NU_PLASMA_FACTOR: f64 = 2. * ELECTRON_CHARGE * ELECTRON_CHARGE / MASS_ELECTRON;

        let inverse_t = -self.0.neg_inverse_t;
        let t = 1. / inverse_t;

        NU_PLASMA_FACTOR
            * theta.sin().powi(2)
            * (inverse_t.besselk(1) + 6. * t * inverse_t.besselk(2))
            / (2. * SPEED_LIGHT * s.powi(2) * inverse_t.besselk(2))
    }

    /// From Heyvaerts+ (2013) equation 43.
    ///
    /// I can get good agreement between this function and Heyvaerts'
    /// isotropic HF computation, *if* T is sufficiently low to make it so
    /// that we can indeed work in the HF approximation.
    fn faraday_v(&self, s: f64, theta: f64) -> f64 {
        const NU_PLASMA_FACTOR: f64 = 2. * ELECTRON_CHARGE * ELECTRON_CHARGE / MASS_ELECTRON;

        let inverse_t = -self.0.neg_inverse_t;

        NU_PLASMA_FACTOR
            * theta.cos()
            * inverse_t.besselk(0)
            / (SPEED_LIGHT * s * inverse_t.besselk(2))
    }
}


#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use ::{Coefficient, DistributionFunction, Stokes, SynchrotronCalculator, gsl};
    use super::ThermalJuettnerDistribution;

    /// Our distributions should integrate to unity in gamma/cos-xi space
    /// when the momentum volume unit of gamma^2*beta is included in the
    /// integrand. The full distribution function in momentum space is our
    /// distribution function divided by (m_e c)^3.
    #[test]
    fn tj_normalization() {
        let mut d = ThermalJuettnerDistribution::new(15.);
        d.normalize();

        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let integral = ws.qagiu(|g| g * (g.powi(2) - 1.).sqrt() * d.calc_f(g, 0.), 1.)
            .tolerance(0., 1e-4)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap();
        assert_approx_eq!(4. * PI * integral, 1., 1e-3);
    }

    /// See the comments in `faraday_q()` — for the approximation to be valid,
    /// *T* must be large and *s* must be very large.
    #[test]
    fn tj_approx_high_freq_rho_q() {
        let mut d = ThermalJuettnerDistribution::new(10.);
        let apx = d.high_freq_approximation();
        let p = apx.compute_dimensionless(Coefficient::Faraday, Stokes::Q, 4e4, 0.4);
        const EXPECTED: f64 = 4.8081e-11;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    #[test]
    fn tj_heyvaerts_high_freq_rho_q() {
        let d = ThermalJuettnerDistribution::new(10.);
        let full = d.full_calculation();
        let p = full.compute_dimensionless(Coefficient::Faraday, Stokes::Q, 4e4, 0.4);
        const EXPECTED: f64 = 4.8081e-11;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    /// Conversely, for Faraday V, we need a low *T* to make the HF
    /// approximation valid.
    #[test]
    fn tj_approx_high_freq_rho_v() {
        let mut d = ThermalJuettnerDistribution::new(0.1);
        let apx = d.high_freq_approximation();
        let p = apx.compute_dimensionless(Coefficient::Faraday, Stokes::V, 40., 0.5);
        const EXPECTED: f64 = 3.064e-4;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    #[test]
    fn tj_heyvaerts_high_freq_rho_v() {
        let d = ThermalJuettnerDistribution::new(0.1);
        let full = d.full_calculation();
        let p = full.compute_dimensionless(Coefficient::Faraday, Stokes::V, 40., 0.5);
        const EXPECTED: f64 = 3.064e-4;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }
}
