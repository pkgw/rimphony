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

use std::f64;

mod gsl;

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


/// The FullSunchrotronCalculator implements the fully detailed
/// double-integral calculation.
#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub struct FullSynchrotronCalculator<D>(D);

#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
enum StokesVSwitch {
    Inactive,
    PositiveLobe,
    NegativeLobe,
}

#[derive(Copy,Clone,Debug,PartialEq)]
struct FullCalculationState<'a, D: 'a> {
    d: &'a D,
    coeff: Coefficient,
    stokes: Stokes,
    s: f64,
    cos_observer_angle: f64,
    sin_observer_angle: f64,
    stokes_v_switch: StokesVSwitch,
}

impl<D: DistributionFunction> SynchrotronCalculator for FullSynchrotronCalculator<D> {
    fn compute_dimensionless(&self, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> f64 {
        FullCalculationState {
            d: &self.0,
            coeff: coeff,
            stokes: stokes,
            s: s,
            cos_observer_angle: theta.cos(),
            sin_observer_angle: theta.sin(),
            stokes_v_switch: StokesVSwitch::Inactive,
        }.compute()
    }
}


impl<'a, D: 'a + DistributionFunction> FullCalculationState<'a, D> {
    fn compute(&mut self) -> f64 {
        const N_MAX: f64 = 30.;

        // This function used to be "n_summation" in symphony's integrate.c.

        let mut ans = 0_f64;

        // Calculate the contributions from the first 30 n's discretely.

        let n_minus = self.s * self.sin_observer_angle.abs();

        self.stokes_v_switch = StokesVSwitch::Inactive;
        let mut gamma_workspace = gsl::IntegrationWorkspace::new(5000);

        for n in ((n_minus + 1.) as i64)..((n_minus + 1. + N_MAX) as i64) {
            ans += self.gamma_integral(&mut gamma_workspace, n as f64);
        }

        // Now integrate the remaining n's, pretending that n can assume any
        // real value. Stokes V is hard to resolve, so if we're computing it,
        // we split the integrals into two lobes and add the results together
        // (see n_integration). We're a bit sloppy about the integration so if
        // it results in a NAN error, just ignore it.

        self.stokes_v_switch = StokesVSwitch::PositiveLobe;

        let n_start = (n_minus + 1. + N_MAX).floor();
        let n_integral_contrib = self.n_integration(n_start).unwrap_or(f64::NAN);

        if !n_integral_contrib.is_nan() {
            ans += n_integral_contrib;
        }

        if self.stokes == Stokes::V {
            self.stokes_v_switch = StokesVSwitch::NegativeLobe;
            let n_integral_contrib = self.n_integration(n_start).unwrap_or(f64::NAN);
            if !n_integral_contrib.is_nan() {
                ans += n_integral_contrib;
            }
        }

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

        ans * match self.coeff {
            Coefficient::Emission => {
                (TWO_PI * ELECTRON_CHARGE).powi(2) /
                    (SPEED_LIGHT * self.cos_observer_angle.abs())
            },
            Coefficient::Absorption => {
                -1. * (TWO_PI * ELECTRON_CHARGE).powi(2) /
                    (2. * MASS_ELECTRON * SPEED_LIGHT * self.cos_observer_angle.abs())
            },
        }
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
        let mut incr_step_factor = 100_f64;
        const DERIV_TOL: f64 = 1e-5;
        const TOLERANCE: f64 = 1e5;

        let mut n_workspace = gsl::IntegrationWorkspace::new(1000);
        let mut gamma_workspace = gsl::IntegrationWorkspace::new(5000);

        // At low harmonic numbers, step conservatively since every n counts.

        if self.s < 10. {
            delta_n = 1.;
            incr_step_factor = 2.;
        }

        while contrib.abs() >= (ans / TOLERANCE).abs() {
            // Evaluate the derivative of the gamma integral with regards to n
            // to figure out the size of the steps we should be taking. Here
            // we've modified Symphony's algorithm to compare `deriv` to
            // `contrib`; otherwise its values are dimensional and can change
            // substantially if we change the units of the integrand.

            let deriv = gsl::deriv_central(|n| self.gamma_integral(&mut gamma_workspace, n), n_start, 1e-8)
                .map(|r| r.value)?;

            if contrib != 0. && (deriv / contrib).abs() < DERIV_TOL {
                delta_n *= incr_step_factor;
            }

            // Compute the next chunk: the integral over all gammas between
            // `n_start` and `n_start + delta_n`.

            contrib = n_workspace.qag(|n| self.gamma_integral(&mut gamma_workspace, n),
                                      n_start, n_start + delta_n)
                .tolerance(0., 1e-3)
                .rule(gsl::IntegrationRule::GaussKonrod31)
                .compute()
                .map(|r| r.value)?;

            ans += contrib;
            n_start += delta_n;
        }

        Ok(ans)
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

        if self.stokes == Stokes::V && self.stokes_v_switch != StokesVSwitch::Inactive {
            if self.stokes_v_switch == StokesVSwitch::PositiveLobe {
                self._gamma_integral_inner(workspace, gamma_peak, gamma_plus_high, n)
            } else {
                self._gamma_integral_inner(workspace, gamma_minus_high, gamma_peak, n)
            }
        } else {
            self._gamma_integral_inner(workspace, gamma_minus_high, gamma_plus_high, n)
        }
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

        gamma * gamma * pol_term * match self.coeff {
            Coefficient::Emission => self.d.calc_f(gamma, cos_xi),
            Coefficient::Absorption => {
                let (dfdg, dfdcx) = self.d.calc_f_derivatives(gamma, cos_xi);
                let dfdcx_factor = (beta * self.cos_observer_angle - cos_xi) / (gamma - 1. / gamma);
                dfdg + dfdcx_factor * dfdcx
            },
        }
    }
}
