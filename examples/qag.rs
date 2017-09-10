extern crate rimphony;

use rimphony::gsl;
use std::f64;

fn inner(i: i32) {
    let mut ws = gsl::IntegrationWorkspace::new(1024);

    let closure = |x: f64| { x.cos().powi(i) };

    println!("Sample: {}", closure(0.));

    let r = ws.qag(closure, 0., 0.5 * f64::consts::PI)
        .tolerance(0., 1e-6)
        .compute();

    println!("Result: {}   Abs err: {}", r.value, r.abserr);
}

fn main() {
    inner(2);
}
