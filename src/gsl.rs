use gsl_sys;
use std::f64;
use std::mem;
use std::os::raw;


#[derive(Debug)]
pub struct IntegrationWorkspace {
    handle: *mut gsl_sys::gsl_integration_workspace,
}

impl IntegrationWorkspace {
    pub fn new(n: usize) -> Self {
        IntegrationWorkspace {
            handle: unsafe { gsl_sys::gsl_integration_workspace_alloc(n) },
        }
    }
}

impl Drop for IntegrationWorkspace {
    fn drop(&mut self) {
        unsafe { gsl_sys::gsl_integration_workspace_free(self.handle) }
    }
}


#[derive(Clone,Copy,Debug,PartialEq)]
pub struct IntegrationResult {
    pub value: f64,
    pub abserr: f64
}


// This enum is unnamed in GSL so I'm not sure how to get bindgen to pick it up.
#[derive(Clone,Copy,Debug,Eq,PartialEq)]
pub enum IntegrationRule {
    GaussKonrod15 = 1,
    GaussKonrod21 = 2,
    GaussKonrod31 = 3,
    GaussKonrod41 = 4,
    GaussKonrod51 = 5,
    GaussKonrod61 = 6,
}


#[derive(Clone,Copy,Debug,Eq,PartialEq)]
pub enum Integrator {
    QAG,
    QAGIU,
}


pub struct IntegrationBuilder<'a, F: 'a> where F: FnMut(f64) -> f64 {
    workspace: &'a mut IntegrationWorkspace,
    function: F,
    kind: Integrator,
    rule: IntegrationRule,
    lower_bound: f64,
    upper_bound: f64,
    epsabs: f64,
    epsrel: f64,
}


extern "C" fn rust_integrand<F>(x: f64, ctxt: *mut raw::c_void) -> f64 where F: FnMut(f64) -> f64 {
    let f: &mut &mut F = unsafe { mem::transmute(ctxt) };
    f(x)
}


impl<'a, F: 'a> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
    pub fn new(ws: &'a mut IntegrationWorkspace, f: F, kind: Integrator,
           lower: f64, upper: f64) -> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
        IntegrationBuilder {
            workspace: ws,
            function: f,
            kind: kind,
            rule: IntegrationRule::GaussKonrod31,
            lower_bound: lower,
            upper_bound: upper,
            epsabs: f64::NAN,
            epsrel: f64::NAN,
        }
    }

    pub fn tolerance(mut self, epsabs: f64, epsrel: f64) -> Self {
        self.epsabs = epsabs;
        self.epsrel = epsrel;
        self
    }

    pub fn rule(mut self, rule: IntegrationRule) -> Self {
        self.rule = rule;
        self
    }

    pub fn compute(mut self) -> IntegrationResult {
        let mut result = f64::NAN;
        let mut abserr = f64::NAN;

        let mut f = gsl_sys::gsl_function_struct {
            function: Some(rust_integrand::<F>),
            params: &mut self.function as *mut _ as *mut raw::c_void,
        };

        let err = match self.kind {
            Integrator::QAG => unsafe {
                gsl_sys::gsl_integration_qag(
                    &mut f,
                    self.lower_bound,
                    self.upper_bound,
                    self.epsabs,
                    self.epsrel,
                    (*self.workspace.handle).size,
                    self.rule as raw::c_int,
                    self.workspace.handle,
                    &mut result,
                    &mut abserr
                )
            },
            Integrator::QAGIU => unsafe {
                gsl_sys::gsl_integration_qagiu(
                    &mut f,
                    self.lower_bound,
                    self.epsabs,
                    self.epsrel,
                    (*self.workspace.handle).size,
                    self.workspace.handle,
                    &mut result,
                    &mut abserr
                )
            },
        };
        // FIXME: don't ignore err

        IntegrationResult { value: result, abserr: abserr }
    }
}


impl IntegrationWorkspace {
    pub fn qagiu<'a, F>(&'a mut self, f: F, lower_bound: f64) -> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
        IntegrationBuilder::new(self, f, Integrator::QAGIU, lower_bound, f64::NAN)
    }

    pub fn qag<'a, F>(&'a mut self, f: F, lower_bound: f64, upper_bound: f64) -> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
        IntegrationBuilder::new(self, f, Integrator::QAG, lower_bound, upper_bound)
    }
}


/// We abuse the IntegrationResult type here since it is exactly what's called
/// for here.
pub fn deriv_central<F>(mut f: F, x: f64, h: f64) -> IntegrationResult where F: FnMut(f64) -> f64 {
    let mut f = gsl_sys::gsl_function_struct {
        function: Some(rust_integrand::<F>),
        params: &mut f as *mut _ as *mut raw::c_void,
    };

    let mut result = f64::NAN;
    let mut abserr = f64::NAN;

    let err = unsafe { gsl_sys::gsl_deriv_central(&mut f, x, h, &mut result, &mut abserr) };
    // FIXME: don't ignore err

    IntegrationResult { value: result, abserr: abserr }
}
