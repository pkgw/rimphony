extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-search=native=/a/lib");
    println!("cargo:rustc-link-lib=dylib=gsl");
    println!("cargo:rustc-link-lib=dylib=gslcblas");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .clang_arg("-I/a/include") // heyoooo
        .whitelisted_type("gsl_.*")
        .whitelisted_function("gsl_.*")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
