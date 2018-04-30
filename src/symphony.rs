// Copyright 2017-2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! Use the "symphony" approach to calculate synchrotron emission and absorption coefficients.

"Symphony" is a synchrotron code written in C that is primarily the work of
Alex Pandya and the group of Charles Gammie at the University of Illinois. The
key publications are [Pandya, Zhang, Chandra, and Gammie (2016;
DOI:10.3847/0004-637X/822/1/34](https://dx.doi.org/10.3847/0004-637X/822/1/34)
and [Leung, Gammie, and Noble
(2011)](https://dx.doi.org/10.1088/0004-637X/737/1/21).

*/

use std::f64;

use gsl;
use leung_bessel;
use slog::Logger;
use super::{Coefficient, DistributionFunction, Stokes};
use super::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI};


#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
enum StokesVSwitch {
    Inactive,
    PositiveLobe,
    NegativeLobe,
}

#[derive(Clone, Debug)]
struct CalculationState<'a, D: 'a> {
    d: &'a D,
    logger: &'a Logger,
    coeff: Coefficient,
    stokes: Stokes,
    s: f64,
    cos_observer_angle: f64,
    sin_observer_angle: f64,
    stokes_v_switch: StokesVSwitch,
}


pub fn compute_dimensionless<'a, D: DistributionFunction>(
    distrib: &'a D, logger: &'a Logger, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64
) -> f64 {
    CalculationState::new(distrib, logger, coeff, stokes, s, theta).compute()
}


impl<'a, D: 'a + DistributionFunction> CalculationState<'a, D> {
    pub fn new(
        distrib: &'a D, logger: &'a Logger, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64
    ) -> Self {
        CalculationState {
            d: distrib,
            logger: logger,
            coeff: coeff,
            stokes: stokes,
            s: s,
            cos_observer_angle: theta.cos(),
            sin_observer_angle: theta.sin(),
            stokes_v_switch: StokesVSwitch::Inactive,
        }
    }

    pub fn compute(&mut self) -> f64 {
        const N_MAX: f64 = 30.;

        trace!(self.logger, "beginning Symphony computation";
               "distrib" => ?self.d,
               "coeff" => ?self.coeff,
               "stokes" => ?self.stokes,
               "s" => self.s,
               "cos_theta" => self.cos_observer_angle,
        );

        // This function used to be "n_summation" in symphony's integrate.c.

        let mut ans = 0_f64;

        // Calculate the contributions from the first 30 n's discretely.

        let n_minus = self.s * self.sin_observer_angle.abs();

        self.stokes_v_switch = StokesVSwitch::Inactive;
        let mut gamma_workspace = gsl::IntegrationWorkspace::new(5000);

        for n in ((n_minus + 1.) as i64)..((n_minus + 1. + N_MAX) as i64) {
            let contrib = self.gamma_integral(&mut gamma_workspace, n as f64);
            ans += contrib;
            trace!(self.logger, "discrete N"; "n" => n, "contrib" => contrib, "ans" => ans);
        }

        // Now integrate the remaining n's, pretending that n can assume any
        // real value. Stokes V is hard to resolve, so if we're computing it,
        // we split the integrals into two lobes and add the results together
        // (see n_integration). We're a bit sloppy about the integration so if
        // it results in a NAN error, just ignore it.

        self.stokes_v_switch = StokesVSwitch::PositiveLobe;

        let n_start = (n_minus + 1. + N_MAX).floor();
        let n_integral_contrib = self.n_integration(n_start).unwrap_or(f64::NAN);
        trace!(self.logger, "N integration 1"; "n_start" => n_start, "contrib" => n_integral_contrib);

        if !n_integral_contrib.is_nan() {
            ans += n_integral_contrib;
        }

        if self.stokes == Stokes::V {
            self.stokes_v_switch = StokesVSwitch::NegativeLobe;
            let n_integral_contrib = self.n_integration(n_start).unwrap_or(f64::NAN);
            trace!(self.logger, "N integration 2"; "n_start" => n_start, "contrib" => n_integral_contrib);
            if !n_integral_contrib.is_nan() {
                ans += n_integral_contrib;
            }
        }

        trace!(self.logger, "answer after Ns"; "ans" => ans);

        // Finally, apply the dimensional constants that don't vary inside any of the
        // integrals.
        //
        // There's a term of `gamma^2 beta m_e^3 c^3` that plays a role in
        // doing a unit conversion in the integrals. We divide out the
        // non-constant parts of this term inside the distribution function
        // and multiply it back in inside the integral. **NOTE** that for
        // absorption coefficients we take the derivative with regards to
        // gamma in between these steps, which is why we don't just remove the
        // superfluous operations. For now.
        //
        // There is also a term of `2 pi` that comes from assuming gyrotropy,
        // i.e. uniformity in electron azimuth angle.
        //
        // There is also a term of 1. / (nu beta |cos theta|) that comes from
        // integrating over the delta function of the resonance condition; see Leung
        // et al. (2011), equation 62.
        //
        // For emission, the physics-based prefactor is `2 pi e^2 nu^2 / c`.
        //
        // For absorption, there's a physics-based prefactor of `-c e^2 / 2 nu`
        // (Pandya equation 12) and one of `2 pi nu / me c^2` (Pandya equation 13).
        //
        // Collecting the terms that aren't gamma or beta:

        let ans = ans * match self.coeff {
            Coefficient::Emission => {
                (TWO_PI * ELECTRON_CHARGE).powi(2) /
                    (SPEED_LIGHT * self.cos_observer_angle.abs())
            },
            Coefficient::Absorption => {
                -1. * (TWO_PI * ELECTRON_CHARGE).powi(2) /
                    (2. * MASS_ELECTRON * SPEED_LIGHT * self.cos_observer_angle.abs())
            },
            Coefficient::Faraday => panic!("should never be reached"),
        };

        trace!(self.logger, "final dimensional Symphony result"; "ans" => ans);
        ans
    }

    /// Compute the contributions to the integration from the region in the
    /// n/gamma plane starting at `n_start` and extending to infinite `n`.
    /// We use an adaptive algorithm that tackles the problem in chunks.
    ///
    /// Symphony had special code for the kappa distribution that identified
    /// the most important `n` values analytically. We've removed the kappa
    /// distribution so that code is gone too.
    fn n_integration(&mut self, mut n_start: f64) -> gsl::GslResult<f64> {
        let mut ans = 0_f64;
        let mut contrib = 0_f64;
        let mut delta_n = 1e5_f64;
        let mut incr_step_factor = 10_f64;
        const DERIV_TOL: f64 = 1e-5;
        const TOLERANCE: f64 = 1e5;

        // Note: if updating numerics, check whether `diagnostic_n_integral`
        // needs to be updated to stay in sync.

        let mut n_workspace = gsl::IntegrationWorkspace::new(1000);
        let mut gamma_workspace = gsl::IntegrationWorkspace::new(5000);

        // At low harmonic numbers, step conservatively since every n counts.

        if self.s < 10. {
            delta_n = 1.;
            incr_step_factor = 2.;
        }

        trace!(self.logger, ". beginning N integration";
               "n_start" => n_start,
               "delta_n" => delta_n,
               "incr_step_factor" => incr_step_factor,
               "DERIV_TOL" => DERIV_TOL,
               "TOLERANCE" => TOLERANCE
        );

        while contrib.abs() >= (ans / TOLERANCE).abs() {
            // Evaluate the derivative of the gamma integral with regards to n
            // to figure out the size of the steps we should be taking. Here
            // we've modified Symphony's algorithm to compare `deriv` to
            // `contrib`; otherwise its values are dimensional and can change
            // substantially if we change the units of the integrand. If the
            // only contributions come from very high n's, both `contrib` and
            // `deriv` are zero, in which case we need to increase our step
            // size.

            let deriv = gsl::deriv_central(|n| self.gamma_integral(&mut gamma_workspace, n), n_start, 1e-8)
                .map(|r| r.value)?;
            trace!(self.logger, ". derivative"; "n_start" => n_start, "deriv" => deriv);

            if deriv == 0. || (contrib != 0. && (deriv / contrib).abs() < DERIV_TOL) {
                delta_n *= incr_step_factor;
                trace!(self.logger, ". increasing delta_n (deriv)"; "delta_n" => delta_n);
            }

            // This is another change from Symphony that vastly speeds up
            // computations for small `s`. It is often the case that the small
            // `delta_n` and `incr_step_factor` used above are far too
            // conservative and make us do way too many integrals. Here we
            // prevent `delta_n` from becoming small compared to `n_start`,
            // which should always be safe.

            if delta_n < n_start / incr_step_factor {
                delta_n *= incr_step_factor;
                trace!(self.logger, ". increasing delta_n (abs)"; "delta_n" => delta_n);
            }

            // Compute the next chunk: the integral over all gammas between
            // `n_start` and `n_start + delta_n`.

            trace!(self.logger, ". integrating gamma_integral"; "lo" => n_start, "hi" => n_start + delta_n);
            contrib = n_workspace.qag(|n| self.gamma_integral(&mut gamma_workspace, n),
                                      n_start, n_start + delta_n)
                .tolerance(0., 1e-3)
                .rule(gsl::IntegrationRule::GaussKonrod31)
                .compute()
                .map(|r| r.value)?;

            ans += contrib;
            trace!(self.logger, ". N contribution"; "contrib" => contrib, "ans" => ans);
            n_start += delta_n;
        }

        Ok(ans)
    }

    /// An internal diagnostic for the N integral.
    fn diagnostic_n_integral(&mut self, n_lo: f64, n_hi: f64) -> gsl::GslResult<f64> {
        let mut n_workspace = gsl::IntegrationWorkspace::new(1000);
        let mut gamma_workspace = gsl::IntegrationWorkspace::new(5000);

        n_workspace.qag(|n| self.gamma_integral(&mut gamma_workspace, n), n_lo, n_hi)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
    }

    /// Integrate across different values of gamma at a fixed single value of
    /// *n*. The equations are such that we can identify where the maximum of
    /// the integral should be and its approximate width.
    fn gamma_integral(&mut self, workspace: &mut gsl::IntegrationWorkspace, n: f64) -> f64 {
        /* Formerly `gamma_integration_result` from Symphony's integrate.c */

        let gamma_minus = (n / self.s -
                           self.cos_observer_angle.abs() *
                           ((n / self.s).powi(2) - self.sin_observer_angle.powi(2)).sqrt()
        ) / self.sin_observer_angle.powi(2);

        let gamma_plus = (n / self.s +
                           self.cos_observer_angle.abs() *
                           ((n / self.s).powi(2) - self.sin_observer_angle.powi(2)).sqrt()
        ) / self.sin_observer_angle.powi(2);

        let gamma_peak = 0.5 * (gamma_plus + gamma_minus);

        const NU_HIGH: f64 = 3e8;
        const NU_LOW: f64 = 1e6;

        let width = if self.s > NU_HIGH {
            1000.
        } else if self.s > NU_LOW {
            10.
        } else {
            1.
        };

        let gamma_minus_high = gamma_peak - (gamma_peak - gamma_minus) / width;
        let gamma_plus_high = gamma_peak - (gamma_peak - gamma_plus) / width;

        //trace!(self.logger, ". . gamma integral";
        //       "n" => n,
        //       "gamma_minus" => gamma_minus,
        //       "gamma_plus" => gamma_plus,
        //       "gamma_peak" => gamma_peak,
        //       "gamma_minus_high" => gamma_minus_high,
        //       "gamma_plus_high" => gamma_plus_high,
        //       "width" => width,
        //);

        let (gamma0, gamma1) = if self.stokes == Stokes::V && self.stokes_v_switch != StokesVSwitch::Inactive {
            if self.stokes_v_switch == StokesVSwitch::PositiveLobe {
                (gamma_peak, gamma_plus_high)
            } else {
                (gamma_minus_high, gamma_peak)
            }
        } else {
            (gamma_minus_high, gamma_plus_high)
        };

        let contrib = self._gamma_integral_inner(workspace, gamma0, gamma1, n);
        trace!(self.logger, ". . gamma integral result";
               "n" => n,
               "gamma0" => gamma0,
               "gamma1" => gamma1,
               "contrib" => contrib
        );
        contrib
    }

    /// An internal diagnostic version of `gamma_integral`.
    fn diagnostic_gamma_integral(&mut self, n: f64) -> f64 {
        let mut ws = gsl::IntegrationWorkspace::new(5000);
        self.gamma_integral(&mut ws, n)
    }

    #[inline]
    fn _gamma_integral_inner(&mut self, workspace: &mut gsl::IntegrationWorkspace,
                             min: f64, max: f64, n: f64) -> f64 {
        // Formerly `gamma_integral` in Symphony's integrate.c
        //
        // In Symphony this function had code to avoid aborts if GSL failed to
        // integrate if the harmonic number was high or the observer angle was
        // low. We properly handle errors so we should never abort. Our policy
        // then becomes more generous: we ignore failures and use zeros
        // instead on *any* failure, not just ones in certain regions of
        // parameter space. In Symphony, failures in those regions would lead
        // to aborts.

        workspace.qag(|g| self.gamma_integrand(g, n), min, max)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(0.)
    }

    /// Calculate the contribution at specific values of *gamma* and *n*.
    fn gamma_integrand(&mut self, gamma: f64, n: f64) -> f64 {
        /* from integrands.c */

        // First, compute the term that depends on which polarization we're
        // talking about. NOTE: Leung et al. (2011) has the wrong sign for
        // Stokes V in its equations, where "wrong" means "not the IEEE/IAU
        // convention". We do it right below.

        let beta = (1. - 1. / (gamma * gamma)).sqrt();
        let cos_xi = (gamma - n / self.s) / (gamma * beta * self.cos_observer_angle);
        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let m = (self.cos_observer_angle - beta * cos_xi) / self.sin_observer_angle;
        let big_n = beta * sin_xi;
        let z = self.s * gamma * beta * self.sin_observer_angle * sin_xi;
        let mj = m * leung_bessel::Jn(n, z);
        let njp = big_n * leung_bessel::Jn_prime(n, z);

        let pol_term = match self.stokes {
            Stokes::I => mj * mj + njp * njp,
            Stokes::Q => mj * mj - njp * njp,
            Stokes::V => 2. * mj * njp,
        };

        // The other bit of the computation depends on on whether we're
        // looking at emission or absorption. See the bottom of self.compute()
        // for explanations of the prefactors; almost all of the ones in the
        // papers can be pulled out of the integral.

        let f_term = match self.coeff {
            Coefficient::Emission => self.d.calc_f(gamma, cos_xi),
            Coefficient::Absorption => {
                let (dfdg, dfdcx) = self.d.calc_f_derivatives(gamma, cos_xi);
                let dfdcx_factor = (beta * self.cos_observer_angle - cos_xi) / (gamma - 1. / gamma);
                dfdg + dfdcx_factor * dfdcx
            },
            Coefficient::Faraday => panic!("should never be reached"),
        };

        let r = gamma * gamma * pol_term * f_term;
        //trace!(self.logger, "gamma_integrand";
        //       "gamma" => gamma,
        //       "n" => n,
        //       "beta" => beta,
        //       "cos_xi" => cos_xi,
        //       "m" => m,
        //       "z" => z,
        //       "mj" => mj,
        //       "njp" => njp,
        //       "pol_term" => pol_term,
        //       "f_term" => f_term,
        //);
        r
    }
}


// Diagnostic support

pub fn diagnostic_n_integral<'a, D: DistributionFunction>(
    distrib: &'a D, logger: &'a Logger, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64,
    n_lo: f64, n_hi: f64
) -> gsl::GslResult<f64> {
    CalculationState::new(distrib, logger, coeff, stokes, s, theta).diagnostic_n_integral(n_lo, n_hi)
}

pub fn diagnostic_gamma_integral<'a, D: DistributionFunction>(
    distrib: &'a D, logger: &'a Logger, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64, n: f64
) -> f64 {
    CalculationState::new(distrib, logger, coeff, stokes, s, theta).diagnostic_gamma_integral(n)
}
