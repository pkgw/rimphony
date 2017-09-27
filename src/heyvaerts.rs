// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! Calculate Faraday coefficients using the Heyvaerts formalism.

This module computes Faraday rotation and conversion coefficients using the
formalism developed by [Heyvaerts et al. (2013;
DOI:10.1093/mnras/stt135)](https://dx.doi.org/10.1093/mnras/stt135).

Our Faraday Q is what Heyvaerts calls "h", and our Faraday V is what Heyvaerts
calls "f".

*/

use gsl;
use special_fun::FloatSpecial;
use std::f64;

use super::{Coefficient, DistributionFunction, Stokes};
use super::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT};

const FOUR_OVER_SQRT_3: f64 = 2.3094010767585034; // 4/sqrt(3)
const INVERSE_C: f64 = 1. / SPEED_LIGHT;
const INVERSE_SQRT_3: f64 = 0.5773502691896257; // 1 / sqrt(3)
const SQRT_8_OVER_3: f64 = 0.9428090415820635; // sqrt(8)/3
const THREE_TWO_THIRDS: f64 = 2.080083823051904; // 3**(2/3)
const TWO_OVER_SQRT_3: f64 = 1.1547005383792517; // 2/sqrt(3)


// Bessel function helpers

/// This is K_{1/3}(x) * L_{1/3}(x)
///
/// Here K is the modified Bessel function of the second kind and L is a
/// similar construction. We use the approximation provided by Heyvaerts.
///
/// TODO: document where the approximation is valid.
///
/// TODO: any advantages if we do the numerics as the difference of two squares?
fn k13l13(x: f64) -> f64 {
    const COEFF: f64 = f64::consts::PI * f64::consts::PI / 3.;
    let plus = x.besseli(1./3.);
    let minus = x.besseli(-1./3.);
    COEFF * (minus - plus) * (minus + plus)
}

/// This is K_{2/3}(x) * L_{2/3}(x)
///
/// Here K is the modified Bessel function of the second kind and L is a
/// similar construction. We use the approximation provided by Heyvaerts.
///
/// TODO: document where the approximation is valid.
///
/// TODO: any advantages if we do the numerics as the difference of two squares?
fn k23l23(x: f64) -> f64 {
    const COEFF: f64 = f64::consts::PI * f64::consts::PI / 3.;
    let plus = x.besseli(2./3.);
    let minus = x.besseli(-2./3.);
    COEFF * (minus - plus) * (minus + plus)
}


// The main algorithm

#[derive(Copy,Clone,Debug,PartialEq)]
struct CalculationState<'a, D: 'a> {
    d: &'a D,
    coeff: Coefficient,
    stokes: Stokes,
    s: f64,
    cos_observer_angle: f64,
    sin_observer_angle: f64,
    sigma0: f64,
    sigma0_sq: f64,

    // Variables set inside integrands:

    sigma: f64,
    pomega: f64,
    x: f64,
    gamma: f64,
    /// mu = cosine pitch angle
    mu: f64,
}


pub fn compute_dimensionless<D: DistributionFunction>(distrib: &D, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> f64 {
    let sigma0 = s * theta.sin();

    CalculationState {
        d: distrib,
        coeff: coeff,
        stokes: stokes,
        s: s,
        cos_observer_angle: theta.cos(),
        sin_observer_angle: theta.sin(),
        sigma0: sigma0,
        sigma0_sq: sigma0.powi(2),
        sigma: f64::NAN,
        pomega: f64::NAN,
        x: f64::NAN,
        gamma: f64::NAN,
        mu: f64::NAN,
    }.compute()
}

impl<'a, D: 'a + DistributionFunction> CalculationState<'a, D> {
    pub fn compute(&mut self) -> f64 {
        // XXX TEMP
        if self.stokes == Stokes::V {
            return 0.;
        }

        let mut ows = gsl::IntegrationWorkspace::new(4096); // outer workspace
        let mut iws = gsl::IntegrationWorkspace::new(4096); // inner workspace

        let mut pomega_left = -3. * self.sigma0;
        let mut pomega_right = 3. * self.sigma0;
        let mut delta_left = pomega_right;
        let mut delta_right = pomega_right;
        let mut keep_going = true;
        const TOL: f64 = 1e-5;
        const DELTA_SCALE_FACTOR: f64 = 5.;

        // The initial computation can work out to 0 if we're in the
        // "low-frequency" regime such that our initial pomegas just don't
        // cover any NR area.

        let mut nr_val = self.nr_outer_integral(&mut ows, &mut iws, pomega_left, pomega_right);
        if nr_val.is_nan() {
            return f64::NAN;
        }

        while keep_going {
            if nr_val != 0. {
                // I think this logic is OK regardless of the sign of
                // rel_deriv, although the motivation is based on a view in
                // which it's negative. We don't take the derivative on the
                // first pass since it is often hard to compute at minimal
                // values of sigma.

                let rel_deriv = gsl::deriv_central(|p| self.nr_outer_integrand(&mut iws, p), pomega_right, 1e-6)
                    .map(|r| r.value)
                    .unwrap_or(f64::NAN);

                if rel_deriv == 0. || (1. / (rel_deriv * delta_right)).abs() > DELTA_SCALE_FACTOR {
                    delta_right *= DELTA_SCALE_FACTOR;
                }
            }

            let contrib = self.nr_outer_integral(&mut ows, &mut iws, pomega_right, pomega_right + delta_right);
            if contrib.is_nan() {
                return f64::NAN;
            }

            if nr_val != 0. {
                keep_going = (contrib / nr_val).abs() > TOL;
            }

            nr_val += contrib;
            pomega_right += delta_right;
        }

        keep_going = true;

        while keep_going {
            let rel_deriv = gsl::deriv_central(|p| self.nr_outer_integrand(&mut iws, p), pomega_left, 1e-6)
                .map(|r| r.value)
                .unwrap_or(f64::NAN);

            if rel_deriv == 0. || (1. / (rel_deriv * delta_left)).abs() > DELTA_SCALE_FACTOR {
                delta_left *= DELTA_SCALE_FACTOR;
            }

            let contrib = self.nr_outer_integral(&mut ows, &mut iws, pomega_left - delta_left, pomega_left);
            if contrib.is_nan() {
                return f64::NAN;
            }

            keep_going = (contrib / nr_val).abs() > TOL;
            nr_val += contrib;
            pomega_left -= delta_left;
        }

        // Now the quasiresonant contribution. Here sigma is the outer
        // integration variable, pomega the inner.

        let mut qr_val = 0.;
        let mut sigma_low = self.sigma0.max(INVERSE_SQRT_3 * self.sigma0.powf(1.5));
        let mut delta_sigma = self.sigma0;
        keep_going = true;

        while keep_going {
            if qr_val != 0. {
                let rel_deriv = gsl::deriv_central(|s| self.qr_outer_integrand(&mut iws, s), sigma_low, 1e-6)
                    .map(|r| r.value)
                    .unwrap_or(f64::NAN);

                if rel_deriv == 0. || (1. / (rel_deriv * delta_sigma)).abs() > DELTA_SCALE_FACTOR {
                    delta_sigma *= DELTA_SCALE_FACTOR;
                }
            }

            let contrib = self.qr_outer_integral(&mut ows, &mut iws, sigma_low, sigma_low + delta_sigma);
            if contrib.is_nan() {
                return f64::NAN;
            }

            if qr_val != 0. {
                keep_going = (contrib / qr_val).abs() > TOL;
            }

            qr_val += contrib;
            sigma_low += delta_sigma;
        }

        // Final scalings:

        2. * ELECTRON_CHARGE.powi(2) * (nr_val + qr_val) /
            (MASS_ELECTRON * (self.s * self.sin_observer_angle).powi(2))
    }

    #[inline]
    fn fill_coord_vars(&mut self, sigma: f64, pomega: f64) {
        self.sigma = sigma;
        self.pomega = pomega;
        self.x = (sigma.powi(2) - pomega.powi(2) - self.sigma0_sq).sqrt();
        self.gamma = (sigma - pomega * self.cos_observer_angle) / (self.sigma0 * self.sin_observer_angle);
        self.mu = (sigma * self.cos_observer_angle - pomega)
            / (self.sigma0 * self.sin_observer_angle * (self.gamma.powi(2) - 1.).sqrt());
    }

    #[inline]
    fn nr_outer_integral(&mut self, ows: &mut gsl::IntegrationWorkspace, iws: &mut gsl::IntegrationWorkspace, p1: f64, p2: f64) -> f64 {
        ows.qag(|p| self.nr_outer_integrand(iws, p), p1, p2)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(f64::NAN)
    }

    fn nr_outer_integrand(&mut self, ws: &mut gsl::IntegrationWorkspace, pomega: f64) -> f64 {
        let sigma_min = (pomega.powi(2) + self.sigma0_sq).sqrt();
        let sigma_max = INVERSE_SQRT_3 * sigma_min.powf(1.5);

        if sigma_max <= sigma_min {
            return 0.;
        }

        let inner = match self.stokes {
            Stokes::Q => |s| {
                self.fill_coord_vars(s, pomega);
                self.h_nr_element()
            },
            _ => panic!("not implemented")
        };

        ws.qag(inner, sigma_min, sigma_max)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(f64::NAN)
    }

    #[inline]
    fn qr_outer_integral(&mut self, ows: &mut gsl::IntegrationWorkspace, iws: &mut gsl::IntegrationWorkspace, s1: f64, s2: f64) -> f64 {
        ows.qag(|s| self.qr_outer_integrand(iws, s), s1, s2)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(f64::NAN)
    }

    fn qr_outer_integrand(&mut self, ws: &mut gsl::IntegrationWorkspace, sigma: f64) -> f64 {
        let pomega_max_phys = (THREE_TWO_THIRDS * sigma.powf(4./3.) - self.sigma0_sq).sqrt();
        let pomega_max_qr = (sigma.powi(2) - self.sigma0_sq).sqrt();
        let pomega_max = pomega_max_phys.min(pomega_max_qr);

        let inner = match self.stokes {
            Stokes::Q => |p| {
                self.fill_coord_vars(sigma, p);
                self.h_qr_element()
            },
            _ => panic!("not implemented")
        };

        ws.qag(inner, -pomega_max, pomega_max)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(f64::NAN)
    }

    /// The quasi-resonant (QR) approximation of the "h" (Faraday Q) coefficient.
    ///
    /// Heyvaerts equation 120 with the standard prefactorization that we have
    /// chosen.
    fn h_qr_element(&self) -> f64 {
        let po_sq = self.pomega.powi(2);
        let smxox = (self.sigma - self.x) / self.x;
        let g = SQRT_8_OVER_3 * (self.sigma - self.x).powf(1.5) / self.x.sqrt();

        let t1 = FOUR_OVER_SQRT_3 * self.x.powi(2) * smxox.powi(2) * k23l23(g);
        let t2 = TWO_OVER_SQRT_3 * po_sq * smxox * k13l13(g);
        let t3 = -f64::consts::PI * (2. * po_sq + self.sigma0_sq) / (po_sq + self.sigma0_sq).sqrt();

        let dfds = self.dfdsigma();

        return INVERSE_C * (t1 + t2 + t3) * dfds;
    }

    /// The nonresonant (NR) approximation of the "h" (Faraday Q) coefficient.
    ///
    /// Heyvaerts equation 119 with the standard prefactorization that we have
    /// chosen.
    fn h_nr_element(&self) -> f64 {
        let s_sq = self.sigma.powi(2);
        let x_sq = self.x.powi(2);
        let ssqmxsq = s_sq - x_sq;
        let a1 = 1. / 8. - 5. / 24. * s_sq / ssqmxsq;
        let a2 = 3. / 128. - 77. / 576. * s_sq / ssqmxsq + 385. / 3456. * (s_sq / ssqmxsq).powi(2);
        let xa1p = -5. / 12. * s_sq * x_sq / ssqmxsq.powi(2);
        let t1 = (6. * a2 - a1.powi(2) + xa1p) / ssqmxsq.sqrt() + a1 * x_sq / ssqmxsq.powf(1.5)
            - x_sq.powi(2) / ssqmxsq.powf(2.5) / 8.;
        let t2 = (6. * a2 - a1.powi(2)) / ssqmxsq.powf(1.5);
        let u1 = 2. * t1 - self.sigma0_sq * t2;

        let dfds = self.dfdsigma();

        return f64::consts::PI * INVERSE_C * u1 * dfds;
    }

    /// The derivative of the distribution function with regards to the
    /// Heyvaerts coordinate sigma.
    fn dfdsigma(&self) -> f64 {
        let (dfdg, dfdcxi) = self.d.calc_f_derivatives(self.gamma, self.mu);

        let g_term = dfdg / (self.sigma0 * self.sin_observer_angle);

        let mu_term = if dfdcxi == 0. {
            0.
        } else {
            // I'm sure this could be simplified but this works so let's just
            // run with it for now. I just banged this out in sympy and checked
            // numerically.
            let q = self.sigma - self.pomega * self.cos_observer_angle;
            let r = self.pomega - self.sigma * self.cos_observer_angle;
            let t = self.sigma0 * self.sin_observer_angle;
            let u = q.powi(2) - t.powi(2);
            let dcxi_dsigma = (q * u * self.cos_observer_angle + u * r + r * t.powi(2))
                / (u.powf(1.5) * q);
            dcxi_dsigma * dfdcxi
        };

        g_term + mu_term
    }
}
