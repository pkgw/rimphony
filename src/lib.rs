extern crate gsl_sys;

use std::f64;

pub mod gsl;

pub const MASS_ELECTRON: f64 = 9.1093826e-28;
pub const SPEED_LIGHT: f64 = 2.99792458e10;
pub const ELECTRON_CHARGE: f64 = 4.80320680e-10;
pub const PLANCKS_CONSTANT: f64 = 6.6260693e-27;
pub const N_MAX: usize = 30;


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

    pub fn finish(self, nu: f64, B: f64, n_e: f64, theta: f64) -> SynchrotronCalculator<Self> {
        SynchrotronCalculator {
            nu: nu,
            magnetic_field: B,
            electron_density: n_e,
            observer_angle: theta,
            d: self,
            nu_c: f64::NAN,
            stokes_v_switch: StokesVSwitch::Inactive,
        }
    }
}


enum StokesVSwitch {
    Inactive,
    PositiveLobe,
    NegativeLobe,
}


pub struct SynchrotronCalculator<D> {
    nu: f64,
    magnetic_field: f64,
    electron_density: f64,
    observer_angle: f64,

    d: D,

    nu_c: f64,
    stokes_v_switch: StokesVSwitch,
}


// Implementation

trait DistributionFunction {
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
        let norm_term = 4. * f64::consts::PI;
        let prefactor = (self.d.p - 1.) / (self.d.gamma_min.powf(1. - self.d.p) - self.d.gamma_max.powf(1. - self.d.p));
        let body = gamma.powf(-self.d.p) * (-gamma / self.d.gamma_cutoff).exp();
        norm_term * prefactor * body
    }
}


impl<D> SynchrotronCalculator<D> where Self: DistributionFunction {
    fn compute(&mut self, coeff: Coefficient) -> f64 {
        let mut ans = 0_f64;

        self.nu_c = ELECTRON_CHARGE * self.magnetic_field / (2. * f64::consts::PI * MASS_ELECTRON * SPEED_LIGHT);

        let n_minus = (self.nu / self.nu_c) * self.observer_angle.sin().abs();

        self.stokes_v_switch = StokesVSwitch::Inactive;

        n_minus
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
