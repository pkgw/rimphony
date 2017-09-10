#[cfg(test)] #[macro_use] extern crate assert_approx_eq;
extern crate gsl_sys;

extern {
    fn my_Bessel_J(n: f64, x: f64) -> f64;
}

pub fn bessel(n: f64, x: f64) -> f64 {
    unsafe { my_Bessel_J(n, x) }
}

#[cfg(test)]
mod tests {
    use super::bessel;

    #[test]
    fn some_values() {
        assert_approx_eq!(bessel(0.,  0.),  1.);
        assert_approx_eq!(bessel(5.,  5.),  0.2611405);
        assert_approx_eq!(bessel(0., 17.), -0.1698543);
    }
}
