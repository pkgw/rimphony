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


/// There's no "U" option because our linear polarization basis is defined
/// such that all U terms are zero.
#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub enum Stokes {
    I,
    Q,
    V
}


#[derive(Copy,Clone,Debug,Eq,Hash,PartialEq)]
pub enum Coefficient {
    Emission(Stokes),
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
    gamma_cutoff: f64,

    norm: Option<f64>,
}


impl PowerLawDistribution {
    pub fn new(p: f64) -> Self {
        PowerLawDistribution {
            p: p,
            gamma_min: 1.,
            gamma_max: 1000.,
            gamma_cutoff: 1e7,
            norm: None,
        }
    }

    fn decache(&mut self) {
        self.norm = None
    }

    pub fn gamma_limits(mut self, gamma_min: f64, gamma_max: f64, gamma_cutoff: f64) -> Self {
        self.gamma_min = gamma_min;
        self.gamma_max = gamma_max;
        self.gamma_cutoff = gamma_cutoff;
        self.decache();
        self
    }

    pub fn finish(self, coeff: Coefficient, nu: f64, b: f64, n_e: f64, theta: f64) -> SynchrotronCalculator<Self> {
        SynchrotronCalculator {
            coeff: coeff,
            nu: nu,
            magnetic_field: b,
            electron_density: n_e,
            observer_angle: theta,
            d: self,
            nu_c: f64::NAN,
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
    magnetic_field: f64,
    electron_density: f64,
    observer_angle: f64,

    d: D,

    nu_c: f64,
    stokes_v_switch: StokesVSwitch,
}


// Implementation

pub trait DistributionFunction {
    fn calc_f(&mut self, gamma: f64) -> f64;

    fn calc_f_for_normalization(&mut self, gamma: f64) -> f64 {
        self.calc_f(gamma)
    }
}


impl DistributionFunction for SynchrotronCalculator<PowerLawDistribution> {
    fn calc_f(&mut self, gamma: f64) -> f64 {
        let norm = if let Some(n) = self.d.norm {
            n
        } else {
            let n = 1. / s_normalize_f(self);
            self.d.norm = Some(n);
            n
        };

        let beta = (1. - 1. / (gamma * gamma)).sqrt();

        let prefactor = self.electron_density * (self.d.p - 1.) /
            (self.d.gamma_min.powf(1. - self.d.p) - self.d.gamma_max.powf(1. - self.d.p));;

        let body = gamma.powf(-self.d.p) * (-gamma / self.d.gamma_cutoff).exp();

        norm * prefactor * body / (SPEED_LIGHT.powi(3) * MASS_ELECTRON.powi(3) * gamma * gamma * beta)
    }

    fn calc_f_for_normalization(&mut self, gamma: f64) -> f64 {
        let norm_term = 4. * PI;
        let prefactor = (self.d.p - 1.) / (self.d.gamma_min.powf(1. - self.d.p) - self.d.gamma_max.powf(1. - self.d.p));
        let body = gamma.powf(-self.d.p) * (-gamma / self.d.gamma_cutoff).exp();
        norm_term * prefactor * body
    }
}


impl<D> SynchrotronCalculator<D> where Self: DistributionFunction {
    pub fn compute(&mut self) -> f64 {
        self.s_n_summation()
    }

    fn s_n_summation(&mut self) -> f64 {
        /* from integrate.c */

        let mut ans = 0_f64;

        self.nu_c = ELECTRON_CHARGE * self.magnetic_field / (TWO_PI * MASS_ELECTRON * SPEED_LIGHT);

        let n_minus = (self.nu / self.nu_c) * self.observer_angle.sin().abs();

        self.stokes_v_switch = StokesVSwitch::Inactive;

        for n in ((n_minus + 1.) as i64)..((n_minus + 1. + N_MAX) as i64) {
            ans += self.s_gamma_integration_result(n as f64);
        }

        self.stokes_v_switch = StokesVSwitch::PositiveLobe;

        let n_integral_contrib = self.s_n_integration(n_minus);

        if !n_integral_contrib.is_nan() {
            ans += n_integral_contrib;
        }

        if self.coeff.stokes() == Stokes::V {
            self.stokes_v_switch = StokesVSwitch::NegativeLobe;
            let n_integral_contrib = self.s_n_integration(n_minus);
            if !n_integral_contrib.is_nan() {
                ans += n_integral_contrib;
            }
        }

        ans
    }

    fn s_n_integration(&mut self, n_minus: f64) -> f64 {
        /* from integrate.c */

        let mut n_start = (N_MAX + n_minus + 1.).floor();
        let mut ans = 0_f64;
        let mut contrib = 0_f64;
        let mut delta_n = 1e5_f64;
        let mut incr_step_factor = 100_f64;
        const DERIV_TOL: f64 = 1e-5;
        const TOLERANCE: f64 = 1e5;

        if self.nu / self.nu_c < 10. {
            delta_n = 1.;
            incr_step_factor = 2.;
        }

        while contrib.abs() >= (ans / TOLERANCE).abs() {
            let deriv = self.s_derivative_of_n(n_start);

            if deriv.abs() < DERIV_TOL {
                delta_n *= incr_step_factor;
            }

            contrib = self.s_n_integral(n_start, n_start + delta_n);
            ans += contrib;
            n_start += delta_n;
        }

        ans
    }

    fn s_derivative_of_n(&mut self, n_start: f64) -> f64 {
        /* from integrate.c */
        gsl::deriv_central(|n| self.s_gamma_integration_result(n), n_start, 1e-8).value
    }

    fn s_n_integral(&mut self, min: f64, max: f64) -> f64 {
        /* from integrate.c */

        let mut ws = gsl::IntegrationWorkspace::new(1000);

        // TODO, maybe: disable GSL errors if nu/nu>c > 1e6 or observer_angle < 0.15

        ws.qag(|n| self.s_gamma_integration_result(n), min, max)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .value
    }

    fn s_gamma_integration_result(&mut self, n: f64) -> f64 {
        /* from integrate.c */

        let gamma_minus = ((n * self.nu_c) / self.nu -
                           self.observer_angle.cos().abs() *
                           ((n * self.nu_c / self.nu).powi(2) - self.observer_angle.sin().powi(2)).sqrt()
        ) / self.observer_angle.sin().powi(2);

        let gamma_plus = ((n * self.nu_c) / self.nu +
                           self.observer_angle.cos().abs() *
                           ((n * self.nu_c / self.nu).powi(2) - self.observer_angle.sin().powi(2)).sqrt()
        ) / self.observer_angle.sin().powi(2);

        let gamma_peak = 0.5 * (gamma_plus + gamma_minus);

        const NU_HIGH: f64 = 3e8;
        const NU_LOW: f64 = 1e6;

        let width = if self.nu / self.nu_c > NU_HIGH {
            1000.
        } else if self.nu / self.nu_c > NU_LOW {
            10.
        } else {
            1.
        };

        let gamma_minus_high = gamma_peak - (gamma_peak - gamma_minus) / width;
        let gamma_plus_high = gamma_peak - (gamma_peak - gamma_plus) / width;

        let result = if self.coeff.stokes() == Stokes::V && self.stokes_v_switch != StokesVSwitch::Inactive {
            let neg_result = self.s_gamma_integral(gamma_minus_high, gamma_peak, n);
            let pos_result = self.s_gamma_integral(gamma_peak, gamma_plus_high, n);

            if self.stokes_v_switch == StokesVSwitch::PositiveLobe {
                pos_result
            } else {
                neg_result
            }
        } else {
            self.s_gamma_integral(gamma_minus_high, gamma_plus_high, n)
        };

        if result.is_nan() {
            0.
        } else {
            result
        }
    }

    fn s_gamma_integral(&mut self, min: f64, max: f64, n: f64) -> f64 {
        let mut ws = gsl::IntegrationWorkspace::new(5000);

        // TODO, maybe: disable GSL errors if nu/nu>c > 1e6 or observer_angle < 0.15

        ws.qag(|g| self.s_gamma_integrand(g, n), min, max)
            .tolerance(0., 1e-3)
            .rule(gsl::IntegrationRule::GaussKonrod31)
            .compute()
            .value
    }

    fn s_polarization_term(&mut self, gamma: f64, n: f64) -> f64 {
        /* from integrands.c */
        let beta = (1. - 1. / (gamma * gamma)).sqrt();
        let cos_xi = (gamma * self.nu - n * self.nu_c) / (gamma * self.nu * beta) * self.observer_angle.sin();
        let m = (self.observer_angle.cos() - beta * cos_xi) / self.observer_angle.sin();
        let n = beta * (1. - cos_xi * cos_xi).sqrt();
        let z = self.nu * gamma * beta * self.observer_angle.sin() * (1. - cos_xi * cos_xi).sqrt() / self.nu_c;
        let k_xx = m * m * leung_bessel::Jn(n, z).powi(2);
        let k_yy = n * n * leung_bessel::Jn_prime(n, z).powi(2);

        match self.coeff.stokes() {
            Stokes::I => k_xx + k_yy,
            Stokes::Q => k_xx - k_yy,
            Stokes::V => -2. * m * n * leung_bessel::Jn(n, z) * leung_bessel::Jn_prime(n, z),
        }
    }

    fn s_numerical_differential_of_f(&mut self, gamma: f64) -> f64 {
        /* from distribution_function_common_routines.c */
        const EPSILON: f64 = 3e-4;
        const GYROPHASE_INDEP: f64 = TWO_PI;
        const _TMP: f64 = MASS_ELECTRON * SPEED_LIGHT;
        const D3P_TO_GAMMA: f64 = _TMP * _TMP * _TMP; // can't powi() for a const
        let prefactor = TWO_PI * self.nu / (MASS_ELECTRON * SPEED_LIGHT * SPEED_LIGHT) *
            GYROPHASE_INDEP * D3P_TO_GAMMA;

        let df = if self.calc_f(gamma + EPSILON).is_nan() {
            (self.calc_f(gamma) - self.calc_f(gamma - EPSILON)) / EPSILON
        } else if self.calc_f(gamma - EPSILON).is_nan() {
            (self.calc_f(gamma + EPSILON) - self.calc_f(gamma)) / EPSILON
        } else {
            (self.calc_f(gamma + EPSILON) - self.calc_f(gamma - EPSILON)) / (2. * EPSILON)
        };

        prefactor * df
    }

    fn s_gamma_integrand(&mut self, gamma: f64, n: f64) -> f64 {
        /* from integrands.c */

        let beta = (1. - 1. / (gamma * gamma)).sqrt();

        let ans = match self.coeff {
            Coefficient::Emission(_) => {
                let func_i = TWO_PI * (ELECTRON_CHARGE * self.nu).powi(2) / SPEED_LIGHT
                    * (MASS_ELECTRON * SPEED_LIGHT).powi(3) * gamma * gamma * beta * TWO_PI
                    * self.calc_f(gamma) * self.s_polarization_term(gamma, n);
                let prefactor = 1. / (self.nu * beta * self.observer_angle.cos().abs());
                prefactor * func_i
            },
            Coefficient::Absorption(_) => {
                let prefactor = -SPEED_LIGHT * ELECTRON_CHARGE * ELECTRON_CHARGE / (2. * self.nu);
                prefactor * gamma * gamma * beta * self.s_numerical_differential_of_f(gamma) *
                    self.s_polarization_term(gamma, n) / (self.nu * beta * self.observer_angle.cos().abs())
            }
        };

        if ans.is_finite() {
            ans
        } else {
            0.
        }
    }
}


// s_ prefix means this is a Symphony transcription
fn s_normalize_f<D: DistributionFunction>(d: &mut D) -> f64 {
    let mut ws = gsl::IntegrationWorkspace::new(1000);

    ws.qagiu(|x| d.calc_f_for_normalization(x), 1.)
        .tolerance(0., 1e-8)
        .compute()
        .value
}
