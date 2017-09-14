# rimphony

This is a reimplementation of the
[Symphony](https://github.com/AFD-Illinois/symphony) synchrotron coefficient
code in the [Rust](https://rust-lang.org/) language.

Testing indicates that the code works well but it is still in development. You
are encouraged to use Symphony if possible, citing
[Pandya et al (2016)](https://dx.doi.org/10.3847/0004-637X/822/1/34) and
[Leung et al (2011)](https://dx.doi.org/10.1088/0004-637X/737/1/21).

I ([PKGW](https://newton.cx/~peter/)) wrote this version to help me understand
Symphony better, to try out Rust as a language for scientific programming, and
because I thought that writing the program in Rust would make it a bit easier
to experiment with some new features, such as pitch-angle dependences in the
electron distribution functions. I feel that that is indeed the case.

## Compiling and running it

You need to have the Rust development tools installed.
[Follow these easy instructions](https://www.rust-lang.org/install.html) if
you don’t.

You also need to have [GSL](https://www.gnu.org/software/gsl/) and its
development files installed. This library is available on most OS’s and
software packaging systems.

Once Rust is available, hopefully all that is needed is to run

```
cargo build
```

and then

```
cargo test
```

to run the test programs. The command

```
cargo run --example one-powerlaw-direct
```

will compile and run a test program.

To actually use this code for an application, you would need to write a
program that depended on the Rimphony “crate” (a Rust term), used the relevant
API calls, and did something with the results. For instance, the example
program `crank-out-pitchypl` computes and prints coefficient values for random
input parameters. I use these as training data for a neural network that can
quickly compute approximate synchrotron coefficients for arbitrary input
parameters.


## Structure

The layout of the files in this repository follows standard Rust practices.

- The `gsl-sys` subdirectory contains a subproject that implements some minimal
  bindings to the [GSL](https://www.gnu.org/software/gsl/) numerics library.
- The `leung-bessel` subdirectory contains a subproject that provides access
  to the fast Bessel function evaluator described in Appendix B of
  [Leung et al (2011)](https://dx.doi.org/10.1088/0004-637X/737/1/21).
- The `src` subdirectory contains the main code.
- The `tests` subdirectory contains some simple tests.
- The `benches` subdirectory contains some minimal benchmarking code.

## Credits

Alex Pandya and Mani Chandra are the main authors of Symphony. Symphony builds
on Harmony, primarily written by Po Kin Leung. This Rust reimplementation is
by Peter K. G. Williams. If using this code in research work, cite the papers
mentioned above.

## License and copyright

This code is copyright Peter Williams and collaboraters, and is licensed under
version 3 of the GPL, like Symphony.
