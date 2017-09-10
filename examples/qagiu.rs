extern crate rimphony;

use rimphony::gsl;

fn main() {
    let mut ws = gsl::IntegrationWorkspace::new(1024);

    let r =  ws.qagiu(|x| { 1. / (x * x) }, 0.5)
        .tolerance(0., 1e-6)
        .compute();

    println!("Result: {}   Abs err: {}", r.value, r.abserr);
}
