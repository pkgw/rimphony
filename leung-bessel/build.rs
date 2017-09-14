// Copyright 2017 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

extern crate gcc;

fn main() {
    gcc::Build::new()
        .flag("-Wall")
        .flag("-I/a/include") // heyoooo, part 2
        .file("src/bessel.c")
        .compile("libleung_bessell.a");

    println!("cargo:rerun-if-changed=src/bessel.c");
}
