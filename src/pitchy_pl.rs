/*! A power-law distribution function with a sine-power pitch angle dependence.

This is an anisotropic power-law distribution. The energy dependence is the
same as the power-law distribution, but there is a dependence on `sin(pitch
angle)^k` for a configurable *k*.

*/

use std::f64;

use gsl;
use super::{TWO_PI, DistributionFunction, FullSynchrotronCalculator};


/// Parameters for a power-law electron distribution with a sin(pitch-angle)
/// dependence. See the module-level documentation for details.
pub struct PitchyPowerLawDistribution {
    p: f64,
    k: f64,
    gamma_min: f64,
    gamma_max: f64,
    inv_gamma_cutoff: f64,
    norm: f64,
}


impl DistributionFunction for PitchyPowerLawDistribution {
    fn calc_f(&self, gamma: f64, cos_xi: f64) -> f64 {
        if gamma < self.gamma_min || gamma > self.gamma_max {
            return 0.;
        }

        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let pa_term = sin_xi.powf(self.k);

        let beta = (1. - 1. / (gamma * gamma)).sqrt();
        let gamma_term = gamma.powf(-self.p) * (-gamma * self.inv_gamma_cutoff).exp();

        self.norm * pa_term * gamma_term / (gamma * gamma * beta)
    }

    fn calc_f_derivatives(&self, gamma: f64, cos_xi: f64) -> (f64, f64) {
        if gamma < self.gamma_min || gamma > self.gamma_max {
            return (0., 0.);
        }

        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let pa_term = sin_xi.powf(self.k);

        let beta = (1. - 1. / (gamma * gamma)).sqrt();
        let gamma_term = gamma.powf(-self.p) * (-gamma * self.inv_gamma_cutoff).exp();

        let f = self.norm * pa_term * gamma_term / (gamma * gamma * beta);

        let dfdg = -f * ((self.p + 1.) / gamma + gamma / (gamma * gamma - 1.) + self.inv_gamma_cutoff);
        let dfdcx = f * self.k * cos_xi / sin_xi;
        (dfdg, dfdcx)
    }
}


impl PitchyPowerLawDistribution {
    /// Create a new set of “pitchy” power-law parameters with the specified
    /// power-law index and sin-pitch-angle power.
    ///
    /// The default gamma limits are a minimum of 1, a maximum of 10^12, and a
    /// cutoff at 10^10.
    pub fn new(p: f64, k: f64) -> Self {
        PitchyPowerLawDistribution {
            p: p,
            k: k,
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

    /// Create a SynchrotronCalculator from this set of parameters. The
    /// calculator will use the full, detailed double integral calculation to
    /// evaluate all coefficients.
    pub fn full_calculation(mut self) -> FullSynchrotronCalculator<Self> {
        // We can compute the pitch-angle normalization factor analytically.

        let pa_integral = gsl::hypergeometric_2F1(0.5, -0.5 * self.k, 1.5, 1.);

        // We probably could do the gamma factor analytically, but we're lazy.

        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let gamma_integral = ws.qag(|g| g.powf(-self.p) * (-g * self.inv_gamma_cutoff).exp(),
                                    self.gamma_min, self.gamma_max)
            .tolerance(0., 1e-8)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap();

        self.norm = 1. / (2. * TWO_PI * pa_integral * gamma_integral);

        // Ready to go.
        FullSynchrotronCalculator(self)
    }
}

#[cfg(test)]
mod test {
    use super::PitchyPowerLawDistribution;
    use ::{Coefficient, PowerLawDistribution, Stokes, SynchrotronCalculator};

    // FIXME: basically copied from benches/powerlaw.rs.

    const SS: &[f64] = &[1e0, 1e1, 1e2, 1e3, 1e4];
    const THETAS: &[f64] = &[0.05, 0.430, 0.810, 1.190, 1.5707];
    const PS: &[f64] = &[1.5, 1.75, 2.5, 3.25, 4.];

    const CHOICES: &[usize] = &[
        1, 4, 2,
        3, 1, 0,
        0, 3, 1,
        2, 0, 4,
        4, 2, 3,
    ];

    const N_CALCS_PER_TEST: usize = 2;

    fn test_k_zero(mut base: usize, coeff: Coefficient, stokes: Stokes) {
        for _ in 0..N_CALCS_PER_TEST {
            let s = SS[CHOICES[base + 0]];
            let theta = THETAS[CHOICES[base + 1]];
            let p = PS[CHOICES[base + 2]];

            let pl = PowerLawDistribution::new(p)
                .full_calculation()
                .compute_dimensionless(coeff, stokes, s, theta);
            let us = PitchyPowerLawDistribution::new(p, 0.)
                .full_calculation()
                .compute_dimensionless(coeff, stokes, s, theta);

            assert_approx_eq!(pl, us);

            base += 3;
        }
    }

    #[test]
    fn k_zero_ji() {
        test_k_zero(0, Coefficient::Emission, Stokes::I);
    }

    #[test]
    fn k_zero_ai() {
        test_k_zero(3, Coefficient::Absorption, Stokes::I);
    }

    #[test]
    fn k_zero_jq() {
        test_k_zero(6, Coefficient::Emission, Stokes::Q);
    }

    #[test]
    fn k_zero_aq() {
        test_k_zero(9, Coefficient::Absorption, Stokes::Q);
    }

    #[test]
    fn k_zero_jv() {
        test_k_zero(0, Coefficient::Emission, Stokes::V);
    }

    #[test]
    fn k_zero_av() {
        test_k_zero(3, Coefficient::Absorption, Stokes::V);
    }
}
