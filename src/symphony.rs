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
            // size. We also set the derivative step size to be a fix fraction
            // of n_start, rather than unconditionally 1e-8, since for
            // challenging integrals (very small pitch angles) we can get n's
            // so big that `n_start + 1e-8 = n_start`.

            let deriv = gsl::deriv_central(|n| self.gamma_integral(&mut gamma_workspace, n),
                                           n_start, 1e-10 * n_start)
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
            trace!(self.logger, ". N contribution";
                   "contrib" => format!("{:.16e}", contrib),
                   "ans" => format!("{:.16e}", ans)
            );
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

        let gamma_peak = 0.5 * (gamma_plus + gamma_minus); // = n / (s sin^2(th))

        // Symphony had a simple heuristic for setting (what we call)
        // `rel_width` based on the value of `s`. But actually, width of the
        // peak around `gamma_peak` depends most significantly on `n`. I have
        // evaluated the integrand for various distribution functions and
        // parameter combinations and found that the following formula gives a
        // good value for `rel_width`. I have tried to make it conservatively
        // wide, but I do find that we need to set `rel_width` to 1 for `s <
        // 1e6` to reproduce Symphony's results. TODO: commit the supporting
        // work used to determine this fitting formula.

        let rel_width = if self.s < 1e6 {
            1.
        } else {
            (-0.27 * n.ln() - 0.1).exp()
        };

        let gamma_minus_high = gamma_peak - (gamma_peak - gamma_minus) * rel_width;
        let gamma_plus_high = gamma_peak - (gamma_peak - gamma_plus) * rel_width;

        //trace!(self.logger, ". . gamma integral";
        //       "n" => n,
        //       "gamma_minus" => gamma_minus,
        //       "gamma_plus" => gamma_plus,
        //       "gamma_peak" => gamma_peak,
        //       "gamma_minus_high" => gamma_minus_high,
        //       "gamma_plus_high" => gamma_plus_high,
        //       "rel_width" => rel_width,
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

        // Formerly `gamma_integral` in Symphony's integrate.c
        //
        // In Symphony this function had code to avoid aborts if GSL failed to
        // integrate if the harmonic number was high or the observer angle was
        // low. We properly handle errors so we should never abort. Our policy
        // then becomes more generous: we ignore failures and use zeros
        // instead on *any* failure, not just ones in certain regions of
        // parameter space. In Symphony, failures in those regions would lead
        // to aborts.

        let contrib = workspace.qag(|g| self.gamma_integrand(g, n), gamma0, gamma1)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap_or(0.);

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

    /// Calculate the contribution at specific values of *gamma* and *n*.
    fn gamma_integrand(&mut self, gamma: f64, n: f64) -> f64 {
        /* from integrands.c */

        // First, compute the term that depends on which polarization we're
        // talking about. NOTE: Leung et al. (2011) has the wrong sign for
        // Stokes V in its equations, where "wrong" means "not the IEEE/IAU
        // convention". We do it right below.

        let beta = (1. - 1. / (gamma * gamma)).sqrt();
        let cos_xi = (self.s * gamma - n) / (self.s * gamma * beta * self.cos_observer_angle);
        let sin_xi = (1. - cos_xi * cos_xi).sqrt();
        let m = (self.cos_observer_angle - beta * cos_xi) / self.sin_observer_angle;
        let big_n = beta * sin_xi;

        // At very low pitch angles and high `s` values, we have to evaluate
        // this integrand at very large `gamma` and `n` (e.g.: theta ~ 0.005,
        // s ~ 1e7, gamma ~ 1e9, n ~ 1e12). The above expression for `sin_xi`
        // loses much of the precision in the input parameter `gamma`, and
        // this noise propagates into `z` and the `n - z` difference that
        // matters in the evaluation of the Bessel functions. The the upshot
        // is that the integrand starts jittering at the ~0.1% level with
        // small changes in `gamma`, which can cause the integrator to freak
        // out and give up.
        //
        // In the computation of `z`, the key value that must be computed
        // carefully is the combination `gamma sin(xi)`. The following
        // equation stabilizes the numerics such that `z` varies smoothly with
        // `gamma` in these challenging situations.

        let beta2_costh2 = (beta * self.cos_observer_angle).powi(2);
        let s_on_r = 2. * n / (self.s * (beta2_costh2 - 1.));
        let r = 1. - 1. / beta2_costh2;
        let gamma_sin_xi = (r * (gamma * (gamma + s_on_r)) - (n * n / (self.s * self.s * beta2_costh2))).sqrt();
        let z = self.s * beta * self.sin_observer_angle * gamma_sin_xi;

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
        //       "gamma" => format!("{:.16e}", gamma),
        //       "n" => format!("{:.16e}", n),
        //       "beta" => format!("{:.16e}", beta),
        //       "cos_xi" => format!("{:.16e}", cos_xi),
        //       "m" => format!("{:.16e}", m),
        //       "z" => format!("{:.16e}", z),
        //       "mj" => format!("{:.16e}", mj),
        //       "njp" => format!("{:.16e}", njp),
        //       "pol_term" => format!("{:.16e}", pol_term),
        //       "f_term" => format!("{:.16e}", f_term),
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

pub fn diagnostic_gamma_integrand<'a, D: DistributionFunction>(
    distrib: &'a D, logger: &'a Logger, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64,
    n: f64, gamma: f64
) -> f64 {
    CalculationState::new(distrib, logger, coeff, stokes, s, theta).gamma_integrand(gamma, n)
}
