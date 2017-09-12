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

extern crate gsl_sys;
extern crate leung_bessel;

use std::f64;

pub mod gsl;

pub use f64::consts::PI;
pub const TWO_PI: f64 = 2. * PI;
pub const MASS_ELECTRON: f64 = 9.1093826e-28;
pub const SPEED_LIGHT: f64 = 2.99792458e10;
pub const ELECTRON_CHARGE: f64 = 4.80320680e-10;
pub const PLANCKS_CONSTANT: f64 = 6.6260693e-27;
pub const N_MAX: f64 = 30.;


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
    Emission(Stokes),

    /// The absorption coefficient, in units of inverse centimeters.
    Absorption(Stokes),
}


impl Coefficient {
    pub fn stokes(&self) -> Stokes {
        match self {
            &Coefficient::Emission(x) => x,
            &Coefficient::Absorption(x) => x,
        }
    }
}


// Distributions

pub struct PowerLawDistribution {
    p: f64,

    gamma_min: f64,
    gamma_max: f64,
    inv_gamma_cutoff: f64,

    norm: f64,
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

    pub fn finish(mut self, coeff: Coefficient, nu: f64, b: f64, n_e: f64, theta: f64) -> SynchrotronCalculator<Self> {
        // Compute the normalization factor.

        let mut ws = gsl::IntegrationWorkspace::new(1000);
        let integral = ws.qag(|g| g.powf(-self.p) * (-g * self.inv_gamma_cutoff).exp(),
                              self.gamma_min, self.gamma_max)
            .tolerance(0., 1e-8)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .map(|r| r.value)
            .unwrap();
        self.norm = n_e / (2. * TWO_PI * integral);

        // Ready to go.

        let nu_c = ELECTRON_CHARGE * b / (TWO_PI * MASS_ELECTRON * SPEED_LIGHT);

        SynchrotronCalculator {
            coeff: coeff,
            nu: nu,
            s: nu / nu_c,
            cos_observer_angle: theta.cos(),
            sin_observer_angle: theta.sin(),
            d: self,
            stokes_v_switch: StokesVSwitch::Inactive,
        }
    }
}


#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
enum StokesVSwitch {
    Inactive,
    PositiveLobe,
    NegativeLobe,
}


pub struct SynchrotronCalculator<D> {
    coeff: Coefficient,
    nu: f64,
    s: f64,
    cos_observer_angle: f64,
    sin_observer_angle: f64,

    d: D,

    stokes_v_switch: StokesVSwitch,
}


// Implementation

pub trait DistributionFunction {
    /// This function gives the modified distribution function `d(n_e) /
    /// (gamma^2 beta d(gamma) d(cos xi) d(phi)])`. In all the cases in
    /// Symphony, isotropy and gyrotropy are assumed, so the `d(cos xi)
    /// d(phi)` become 1/4pi. The funny scalings are due to me being stuck on
    /// the numerical derivative step. To be revisited.
    fn calc_f(&mut self, gamma: f64) -> f64;
}


impl DistributionFunction for SynchrotronCalculator<PowerLawDistribution> {
    fn calc_f(&mut self, gamma: f64) -> f64 {
        if gamma < self.d.gamma_min || gamma > self.d.gamma_max {
            0.
        } else {
            let beta = (1. - 1. / (gamma * gamma)).sqrt();

            self.d.norm * gamma.powf(-self.d.p) * (-gamma * self.d.inv_gamma_cutoff).exp()
                / (gamma * gamma * beta)
        }
    }
}


impl<D> SynchrotronCalculator<D> where Self: DistributionFunction {
    pub fn compute(&mut self) -> f64 {
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

        if self.coeff.stokes() == Stokes::V {
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
            Coefficient::Emission(_) => {
                (TWO_PI * ELECTRON_CHARGE).powi(2) * self.nu /
                    (SPEED_LIGHT * self.cos_observer_angle.abs())
            },
            Coefficient::Absorption(_) => {
                -1. * (TWO_PI * ELECTRON_CHARGE).powi(2) /
                    (2. * MASS_ELECTRON * SPEED_LIGHT * self.nu * self.cos_observer_angle.abs())
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
            // Evaluate the derivative of the gamma integral with regards to n to figure out
            // the size of the steps we should be taking.

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

        if self.coeff.stokes() == Stokes::V && self.stokes_v_switch != StokesVSwitch::Inactive {
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

        let pol_term = match self.coeff.stokes() {
            Stokes::I => mj * mj + njp * njp,
            Stokes::Q => mj * mj - njp * njp,
            Stokes::V => 2. * mj * njp,
        };

        // The other bit of the computation depends on on whether we're
        // looking at emission or absorption. See the bottom of self.compute()
        // for explanations of the prefactors; almost all of the ones in the
        // papers can be pulled out of the integral.

        let ans = match self.coeff {
            Coefficient::Emission(_) => gamma * gamma * self.calc_f(gamma) * pol_term,
            Coefficient::Absorption(_) => gamma * gamma * self.numerical_differential_of_f(gamma) * pol_term,
        };

        if ans.is_finite() {
            ans
        } else {
            0.
        }
    }

    fn numerical_differential_of_f(&mut self, gamma: f64) -> f64 {
        /* from distribution_function_common_routines.c */
        const EPSILON: f64 = 3e-4;

        let f_plus = self.calc_f(gamma + EPSILON);

        if f_plus.is_nan() {
            (self.calc_f(gamma) - self.calc_f(gamma - EPSILON)) / EPSILON
        } else {
            let f_minus = self.calc_f(gamma - EPSILON);

            if f_minus.is_nan() {
                (f_plus - self.calc_f(gamma)) / EPSILON
            } else {
                (f_plus - f_minus) / (2. * EPSILON)
            }
        }
    }
}
