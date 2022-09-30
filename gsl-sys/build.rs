// Copyright 2017-2019 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

extern crate bindgen;
extern crate pkg_config;

use std::env;
use std::path::PathBuf;

fn main() {
    let gsl = pkg_config::Config::new()
        .atleast_version("1.0")
        .probe("gsl")
        .unwrap();

    let mut builder = bindgen::Builder::default().header("src/wrapper.h");

    for ref path in &gsl.include_paths {
        builder = builder.clang_arg(format!("-I{}", path.display()));
    }

    let bindings = builder
        .allowlist_type("gsl_.*")
        .allowlist_function("gsl_.*")
        .blocklist_function("gsl_coerce_long_double") // uses u128 which is not stable
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
