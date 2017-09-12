/*! The power-law distribution function.

The electrons are isotropic. The distribution is zero outside of the bounds
`gamma_min` and `gamma_max`. The power-law index is `p`, such that `dN/dgamma
~ gamma^(-p)`. An exponential cutoff of the form `exp(-gamma/gamma_cutoff)` is
multiplied in. This cutoff shoud be smaller than `gamma_max` to prevent the
integrators from having problems with the hard cutoff at `gamma_max`.

*/

use std::f64;

use gsl;
use super::{TWO_PI, DistributionFunction, FullSynchrotronCalculator};


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
    pub fn new(p: f64) -> Self {
        PowerLawDistribution {
            p: p,
            gamma_min: 1.,
            gamma_max: 1e12,
            inv_gamma_cutoff: 1e-10,
            norm: f64::NAN,
        }
    }

    pub fn gamma_limits(mut self, gamma_min: f64, gamma_max: f64, gamma_cutoff: f64) -> Self {
        self.gamma_min = gamma_min;
        self.gamma_max = gamma_max;
        self.inv_gamma_cutoff = 1. / gamma_cutoff;
        self
    }

    pub fn full_calculation(mut self) -> FullSynchrotronCalculator<Self> {
        // Compute the normalization factor.

        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let integral = ws.qag(|g| g.powf(-self.p) * (-g * self.inv_gamma_cutoff).exp(),
                              self.gamma_min, self.gamma_max)
            .tolerance(0., 1e-8)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap();
        self.norm = 1. / (2. * TWO_PI * integral);

        // Ready to go.
        FullSynchrotronCalculator(self)
    }
}
