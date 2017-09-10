#[cfg(test)] #[macro_use] extern crate assert_approx_eq;
extern crate gsl_sys;

extern {
    fn my_Bessel_J(n: f64, x: f64) -> f64;
    fn my_Bessel_dJ(n: f64, x: f64) -> f64;
}

#[allow(non_snake_case)]
pub fn Jn(n: f64, x: f64) -> f64 {
    unsafe { my_Bessel_J(n, x) }
}

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
