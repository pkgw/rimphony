// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

/*!
This crate provides a functions that can compute Bessel functions of large *n*
very quickly.

In particular, it can computes Bessel functions of the first kind with the
order and the argument both real and nonnegative. To learn more about Bessel
functions, try
[Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J.CE.B1)
or [Wolfram
MathWorld](http://mathworld.wolfram.com/BesselFunctionoftheFirstKind.html).

The implementation is by Po Kin Leung and is described in Leung, Gammie, and
Noble (2011,
DOI:[10.1088/0004-637X/737/1/21](https://dx.doi.org/10.1088/0004-637X/737/1/21))
and was made public in their C program “harmony”. See that paper for the gory
details, motivation, and references to the various Bessel function
approximations that are used in different conditions. The implementation used
here is taken from the successor program,
[symphony](https://github.com/AFD-Illinois/symphony). Specifically, it is
based on the version in Git commit `41a42e545b003e7df71ad286b1d4bba2fa329702`.
It has been tidied up by Peter K. G. Williams.

If the order of the Bessel function is less than 30, it must be an integer,
and the GSL implementation is used to compute its value. If the argument is
not integral, NaN is returned.

*/

#[cfg(test)] #[macro_use] extern crate assert_approx_eq;
extern crate gsl_sys;
extern crate libc;

extern {
    fn my_Bessel_J(n: libc::c_double, x: libc::c_double) -> libc::c_double;
    fn my_Bessel_dJ(n: libc::c_double, x: libc::c_double) -> libc::c_double;
}


/// Compute the Bessel function of the first kind, `J_n`.
///
/// The order is *n* and the argument is *x*. See the crate-level
/// documentation for details. Both *n* and *x* must be nonnegative. If *n* is
/// less than 30, it must be an integer, or NaN will be returned.
///
/// To learn more about Bessel functions, try
/// [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J.CE.B1)
/// or [Wolfram
/// MathWorld](http://mathworld.wolfram.com/BesselFunctionoftheFirstKind.html).
#[allow(non_snake_case)]
pub fn Jn(n: f64, x: f64) -> f64 {
    unsafe { my_Bessel_J(n, x) }
}


/// Compute the derivative of the Bessel function of the first kind, `J'_n`.
///
/// The order is *n* and the argument is *x*, and the derivative is with
/// regards to the argument. See the crate-level documentation for details.
/// Both *n* and *x* must be nonnegative. If *n* is less than 30, it must be
/// an integer, or NaN will be returned.
///
/// To learn more about Bessel functions, try
/// [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J.CE.B1)
/// or [Wolfram
/// MathWorld](http://mathworld.wolfram.com/BesselFunctionoftheFirstKind.html).
#[allow(non_snake_case)]
pub fn Jn_prime(n: f64, x: f64) -> f64 {
    unsafe { my_Bessel_dJ(n, x) }
}

#[cfg(test)]
mod tests {
    use super::Jn;

    #[test]
    fn some_values() {
        assert_approx_eq!(Jn(0.,  0.),  1.);
        assert_approx_eq!(Jn(5.,  5.),  0.2611405);
        assert_approx_eq!(Jn(0., 17.), -0.1698543);
    }
}
