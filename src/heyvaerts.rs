// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*! Calculate Faraday coefficients using the Heyvaerts formalism.

This module computes Faraday rotation and conversion coefficients using the
formalism developed by [Heyvaerts et al. (2013;
DOI:10.1093/mnras/stt135)](https://dx.doi.org/10.1093/mnras/stt135).

*/

use std::f64;

use gsl;
use leung_bessel;
use super::{Coefficient, DistributionFunction, Stokes};
use super::{ELECTRON_CHARGE, MASS_ELECTRON, SPEED_LIGHT, TWO_PI};


#[derive(Copy,Clone,Debug,PartialEq)]
pub struct CalculationState<'a, D: 'a> {
    d: &'a D,
    coeff: Coefficient,
    stokes: Stokes,
    s: f64,
    cos_observer_angle: f64,
    sin_observer_angle: f64,
}


impl<'a, D: 'a + DistributionFunction> CalculationState<'a, D> {
    pub fn new(distrib: &'a D, coeff: Coefficient, stokes: Stokes, s: f64, theta: f64) -> Self {
        CalculationState {
            d: distrib,
            coeff: coeff,
            stokes: stokes,
            s: s,
            cos_observer_angle: theta.cos(),
            sin_observer_angle: theta.sin(),
        }
    }

    pub fn compute(&mut self) -> f64 {
        0.
    }
}
