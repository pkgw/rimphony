// Copyright 2017-2018 Peter Williams <peter@newton.cx> and collaborators
// Licensed under the GPL version 3.

//! A tiny helper for testing convenience.

extern crate rand;
#[macro_use] extern crate slog;
extern crate slog_async;
extern crate slog_term;

use slog::Drain;

/// Create a simple `slog` logger for use in test programs.
///
/// It logs to the terminal using default parameters, as per the `slog` basic
/// example. This just saves us ~8 lines of boilerplate in all of our
/// test/demo programs.
pub fn default_log() -> slog::Logger {
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::FullFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain)
        .overflow_strategy(slog_async::OverflowStrategy::Block)
        .build().fuse();
    slog::Logger::root(drain, o!())
}


/// A simple utility for sampling random numbers.
///
/// The distribution can be uniform or log-uniform.
pub struct Sampler {
    is_log: bool,
    low: f64,
    range: f64
}

impl Sampler {
    /// Create a new Sampler.
    pub fn new(is_log: bool, mut low: f64, mut high: f64) -> Self {
        if low > high {
            let tmp = high;
            high = low;
            low = tmp;
        }

        if is_log {
            low = low.ln();
            high = high.ln();
        }

        Sampler { is_log: is_log, low: low, range: high - low }
    }

    /// Sample a number from the distribution.
    pub fn get(&self) -> f64 {
        let n = self.low + rand::random::<f64>() * self.range;

        if self.is_log {
            n.exp()
        } else {
            n
        }
    }
}
