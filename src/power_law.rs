// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! The power-law distribution function.

The electrons are isotropic. The distribution is zero outside of the bounds
`gamma_min` and `gamma_max`. The power-law index is `p`, such that `dN/dgamma
~ gamma^(-p)`. (This means that the shape of the distribution function in
momentum space is *not* a simple power law, because of the nonlinear mapping
between gamma and momentum.) An exponential cutoff of the form
`exp(-gamma/gamma_cutoff)` is multiplied in. This cutoff shoud be smaller than
`gamma_max` to prevent the integrators from having problems with the hard
cutoff at `gamma_max`.

*/

use std::f64;

use gsl;
use super::{TWO_PI, Coefficient, DistributionFunction, FullSynchrotronCalculator, Stokes, SynchrotronCalculator};


/// Parameters for a power-law electron distribution. See the module-level
/// documentation for details.
#[derive(Copy,Clone,Debug,PartialEq)]
pub struct PowerLawDistribution {
    p: f64,
    gamma_min: f64,
    gamma_max: f64,
    inv_gamma_cutoff: f64,
    norm: f64,
}


impl DistributionFunction for PowerLawDistribution {
    fn calc_f(&self, gamma: f64, _cos_xi: f64) -> f64 {
        if gamma < self.gamma_min || gamma > self.gamma_max {
            0.
        } else {
            let beta = (1. - 1. / (gamma * gamma)).sqrt();

            self.norm * gamma.powf(-self.p) * (-gamma * self.inv_gamma_cutoff).exp()
                / (gamma * gamma * beta)
        }
    }

    fn calc_f_derivatives(&self, gamma: f64, _cos_xi: f64) -> (f64, f64) {
        if gamma < self.gamma_min || gamma > self.gamma_max {
            return (0., 0.);
        }

        let p_plus_1 = self.p + 1.;
        let g2_minus_1 = gamma * gamma - 1.;
        let dfdg = -self.norm * gamma.powf(-p_plus_1) / g2_minus_1.sqrt() *
            (-gamma * self.inv_gamma_cutoff).exp() *
            (p_plus_1 / gamma + gamma / g2_minus_1 + self.inv_gamma_cutoff);
        let dfdcx = 0.;

        (dfdg, dfdcx)
    }
}


impl PowerLawDistribution {
    /// Create a new set of power-law parameters with the specified power-law
    /// index.
    ///
    /// The default gamma limits are a minimum of 1, a maximum of 10^12, and a
    /// cutoff at 10^10.
    pub fn new(p: f64) -> Self {
        PowerLawDistribution {
            p: p,
            gamma_min: 1.,
            gamma_max: 1e12,
            inv_gamma_cutoff: 1e-10,
            norm: f64::NAN,
        }
    }

    /// Alter the gamma integration limits of this distribution.
    pub fn gamma_limits(mut self, gamma_min: f64, gamma_max: f64, gamma_cutoff: f64) -> Self {
        self.gamma_min = gamma_min;
        self.gamma_max = gamma_max;
        self.inv_gamma_cutoff = 1. / gamma_cutoff;
        self
    }

    /// Compute the normalization factor. The definition of our distribution
    /// function is such that the integral
    fn normalize(&mut self) {
        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let integral = ws.qag(|g| g.powf(-self.p) * (-g * self.inv_gamma_cutoff).exp(),
                              self.gamma_min, self.gamma_max)
            .tolerance(0., 1e-8)
            .rule(gsl::IntegrationRule::GaussKonrod31)
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
pub struct HighFrequencyApproximation<'a>(&'a PowerLawDistribution);

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
    /// From Huang & Shcherbakov (2011) equation 51.
    ///
    /// To reproduce Figure 6 of Huang & Shcherbakov, use n = 2.5, s = 1e4,
    /// theta = pi/4, omega_p such that n = 1, and omega = 2 pi. For gamma_min
    /// = 10, I get 1.8e-9, which looks about right.
    fn faraday_q(&self, s: f64, theta: f64) -> f64 {
        0.0085 *
            2. / (self.0.p - 2.) *
            ((s / (theta.sin() * self.0.gamma_min.powi(2))).powf((self.0.p - 2.) / 2.) - 1.) *
            (self.0.p - 1.) / self.0.gamma_min.powf(1. - self.0.p) *
            (theta.sin() / s).powf((self.0.p + 2.) / 2.)
    }

    /// From Huang & Shcherbakov (2011) equation 51.
    ///
    /// This one is not so useful for reproducing their Figure 6 since it's
    /// nigh impossible to read accurate quantitative values for rho_V off of
    /// their plot.
    fn faraday_v(&self, s: f64, theta: f64) -> f64 {
        0.017
            * (self.0.gamma_min.ln() * (self.0.p - 1.))
            / ((self.0.p + 1.) * self.0.gamma_min.powi(2))
            / s
            * theta.sin()
    }
}


#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use ::{Coefficient, DistributionFunction, Stokes, SynchrotronCalculator, gsl};
    use super::PowerLawDistribution;

    /// Our distributions should integrated to unity in gamma/cos-xi space
    /// when the momentum volume unit of gamma^2*beta is included in the
    /// integrand. The full distribution function in momentum space is our
    /// distribution function divided by (m_e c)^3.
    #[test]
    fn test_normalization() {
        let mut d = PowerLawDistribution::new(2.5).gamma_limits(10., 1e12, 1e10);
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

    #[test]
    fn test_hs11_high_freq_rho_q() {
        let mut d = PowerLawDistribution::new(2.5).gamma_limits(10., 1e12, 1e10);
        let apx = d.high_freq_approximation();
        let p = apx.compute_dimensionless(Coefficient::Faraday, Stokes::Q, 1e4, 0.25 * PI);
        const EXPECTED: f64 = 1.81e-9;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    #[test]
    fn test_heyvaerts_high_freq_rho_q() {
        let d = PowerLawDistribution::new(2.5).gamma_limits(10., 1e12, 1e10);
        let full = d.full_calculation();
        let p = full.compute_dimensionless(Coefficient::Faraday, Stokes::Q, 1e4, 0.25 * PI);
        const EXPECTED: f64 = 1.89e-9;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    /// Note that this disagrees with the full integration somewhat
    /// substantially, but as best I can convince myself, the full integration
    /// is just more correct here. As per the left panel of HS11 Figure 6, the
    /// approximation does not do a very good job for gamma_min <~ 100 (and
    /// possible even for larger values of gamma_min, but you just can't tell
    /// from the plot.
    #[test]
    fn test_hs11_high_freq_rho_v() {
        let mut d = PowerLawDistribution::new(2.5).gamma_limits(10., 1e12, 1e10);
        let apx = d.high_freq_approximation();
        let p = apx.compute_dimensionless(Coefficient::Faraday, Stokes::V, 1e4, 0.25 * PI);
        const EXPECTED: f64 = 1.19e-8;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }

    #[test]
    fn test_heyvaerts_high_freq_rho_v() {
        let d = PowerLawDistribution::new(2.5).gamma_limits(10., 1e12, 1e10);
        let full = d.full_calculation();
        let p = full.compute_dimensionless(Coefficient::Faraday, Stokes::V, 1e4, 0.25 * PI);
        const EXPECTED: f64 = 5.28e-8;
        assert_approx_eq!(p, EXPECTED, 0.01 * EXPECTED);
    }
}
