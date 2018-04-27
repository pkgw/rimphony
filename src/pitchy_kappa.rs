// Copyright 2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! A kappa distribution function with a sine-power pitch angle dependence.

We use the relativistic kappa distribution as defined in Pandya et al 2016,
which cites Xiao et al 2006 as its reference. Our γ²β denominator cancels out
a term like this in the numerator (as in the thermal Jüttner distribution).

This is an anisotropic kappa distribution. The energy dependence is the
same as a kappa distribution, but there is a dependence on `sin(pitch
angle)^k` for a configurable *k*.

XXX THIS NEEDS MORE TESTING.

*/

use std::f64;

use gsl;
use super::{TWO_PI, DistributionFunction, FullSynchrotronCalculator};


/// Parameters for a kappa electron distribution with a sin(pitch-angle)
/// dependence. See the module-level documentation for details.
#[derive(Copy,Clone,Debug,PartialEq)]
pub struct PitchyKappaDistribution {
    kappa: f64,
    width: f64,
    inv_kappa_width: f64,
    k: f64,
    inv_gamma_cutoff: f64,
    norm: f64,
}


impl DistributionFunction for PitchyKappaDistribution {
    fn calc_f(&self, gamma: f64, cos_xi: f64) -> f64 {
        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let pa_term = sin_xi.powf(self.k);
        let gamma_term =
            (1. + (gamma - 1.) * self.inv_kappa_width).powf(-(self.kappa + 1.)) *
            (-gamma * self.inv_gamma_cutoff).exp();

        self.norm * pa_term * gamma_term
    }

    fn calc_f_derivatives(&self, gamma: f64, cos_xi: f64) -> (f64, f64) {
        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let pa_term = sin_xi.powf(self.k);
        let gamma_term =
            (1. + (gamma - 1.) * self.inv_kappa_width).powf(-(self.kappa + 1.)) *
            (-gamma * self.inv_gamma_cutoff).exp();

        let f = self.norm * pa_term * gamma_term;

        let dfdg = -f * ((self.kappa + 1.) / (self.kappa * self.width + gamma - 1.) + self.inv_gamma_cutoff);
        let dfdcx = -f * self.k * cos_xi / (sin_xi * sin_xi);
        (dfdg, dfdcx)
    }
}


impl PitchyKappaDistribution {
    /// Create a new set of “pitchy” kappa parameters with the specified
    /// kappa, width, and sin-pitch-angle power.
    ///
    /// The default gamma cutoff is 10^10.
    pub fn new(kappa: f64, width: f64, k: f64) -> Self {
        PitchyKappaDistribution {
            kappa: kappa,
            width: width,
            inv_kappa_width: 1. / (kappa * width),
            k: k,
            inv_gamma_cutoff: 1e-10,
            norm: f64::NAN,
        }
    }

    /// Alter the gamma cutoff of this distribution.
    pub fn gamma_cutoff(mut self, gamma_cutoff: f64) -> Self {
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

        let gamma_integral = {
            let fg = |g: f64| {
                g * (g.powi(2) - 1.).sqrt() *
                    (1. + (g - 1.) * self.inv_kappa_width).powf(-(self.kappa + 1.)) *
                    (-g * self.inv_gamma_cutoff).exp()
            };

            ws.qagiu(fg, 1.)
                .tolerance(0., 1e-8)
                .rule(gsl::IntegrationRule::GaussKonrod31)
                .compute()
                .map(|r| r.value)
                .unwrap()
        };

        self.norm = 1. / (2. * TWO_PI * pa_integral * gamma_integral);

        // Ready to go.
        FullSynchrotronCalculator(self)
    }
}

#[cfg(test)]
mod tests {
    use std::f64;

    use super::PitchyKappaDistribution;
    use ::DistributionFunction;

    #[test]
    fn test_derivatives() {
        use rand;
        const EPS: f64 = 1e-6;
        const TOL: f64 = 1e-4;

        for _ in 0..100 {
            let kappa = 1.5 + 3. * rand::random::<f64>();
            let width = (1. + 2. * rand::random::<f64>()).exp();
            let k = 0. + 3. * rand::random::<f64>();
            let gamma = 1.1 + 1e3 * rand::random::<f64>();
            let cos_xi = 0.01 + 0.98 * rand::random::<f64>();

            let mut ppd = PitchyKappaDistribution::new(kappa, width, k);
            ppd.norm = 1.; // fake this

            let (analytic_dfdg, analytic_dfdcx) = ppd.calc_f_derivatives(gamma, cos_xi);

            let f0 = ppd.calc_f(gamma, cos_xi);
            let numeric_dfdg = (ppd.calc_f(gamma + EPS, cos_xi) - f0) / EPS;
            let numeric_dfdcx = (ppd.calc_f(gamma, cos_xi + EPS) - f0) / EPS;

            // The logical nots here let us catch NaNs.

            if !(((analytic_dfdg - numeric_dfdg) / numeric_dfdg).abs() < TOL) {
                panic!("numerical gamma derivative failed: kappa={:.16e} width={:.16e} \
                        k={:.16e} gamma={:.16e} cosxi={:.16e} analytic={:.16e} \
                        numeric={:.16e}", kappa, width, k, gamma, cos_xi, analytic_dfdg,
                       numeric_dfdg);
            }

            if !(((analytic_dfdcx - numeric_dfdcx) / numeric_dfdcx).abs() < TOL) {
                panic!("numerical cosxi derivative failed: kappa={:.16e} width={:.16e} \
                        k={:.16e} gamma={:.16e} cosxi={:.16e} analytic={:.16e} \
                        numeric={:.16e}", kappa, width, k, gamma, cos_xi, analytic_dfdcx,
                       numeric_dfdcx);
            }
        }
    }
}
