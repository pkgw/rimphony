# rimphony

Rimphony is a tool to calculate the eight coefficients needed to do numerical
radiative transfer of polarized synchrotron emission. Unlike some comparable
tools, it supports particle distributions that are anisotropic in pitch angle,
and can calculate Faraday rotation and conversion coefficients. And, of
course, it’s fully open source!

Rimphony began life as an experimental reimplementation of the
[Symphony](https://github.com/AFD-Illinois/symphony) synchrotron coefficient
code in the [Rust](https://rust-lang.org/) language. Symphony comes out of
Charles Gammie's group and is primarily the work of Alex Pandya. Its
algorithms are described in
[Pandya et al (2016)](https://dx.doi.org/10.3847/0004-637X/822/1/34) and
[Leung et al (2011)](https://dx.doi.org/10.1088/0004-637X/737/1/21).

Rimphony has adapted the computational framework of Symphony and adds new
routines to calculate Faraday rotation and conversion coefficients using the
formalism developed by
[Heyvaerts et al (2013; DOI:10.1093/mnras/stt135)](https://dx.doi.org/10.1093/mnras/stt135).
Some of the Symphony implementation details have also been altered to increase
speed and reliability in challenging corners of parameter space. While
rimphony is not yet described in a refereed publication, based on the results
of its built-in validation test suite and benchmarking, I
([PKGW](https://newton.cx/~peter/)) feel comfortable recommending it for use
in research projects.

As a side note, my overall takeaway is that Rust is an *excellent* language
for this kind of project. In my opinion, the code of Rimphony is substantially
more readable and extensible than that of Symphony, and the infrastructure of
the Rust language make compilation, testing, benchmarking, documentation, and
distribution all straightforward.


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
parameters — the [neurosynchro](https://github.com/pkgw/neurosynchro/)
project.


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
- The `examples` subdirectory contains a few example programs that exercise
  the library.


## Credits

Alex Pandya and Mani Chandra are the main authors of Symphony. Symphony builds
on Harmony, primarily written by Po Kin Leung. This Rust reimplementation is
by Peter K. G. Williams. If using this code in research work, cite the papers
mentioned above.


## License and copyright

This code is copyright Peter Williams and collaborators, and is licensed under
version 3 of the GPL, like Symphony.
