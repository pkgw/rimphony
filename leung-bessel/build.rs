// Copyright 2017-2019 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

extern crate cc;

fn main() {
    cc::Build::new()
        .flag("-Wall")
        .flag("-I/a/include") // heyoooo, part 2
        .file("src/bessel.c")
        .compile("libleung_bessell.a");

    println!("cargo:rerun-if-changed=src/bessel.c");
}
