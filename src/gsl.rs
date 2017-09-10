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


#[derive(Clone,Copy,Debug,Eq,PartialEq)]
enum Integrator {
    QAGIU,
}


pub struct IntegrationBuilder<'a, F: 'a> where F: FnMut(f64) -> f64 {
    workspace: &'a mut IntegrationWorkspace,
    function: F,
    kind: Integrator,
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
    fn new(ws: &'a mut IntegrationWorkspace, f: F, kind: Integrator,
           lower: f64, upper: f64) -> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
        IntegrationBuilder {
            workspace: ws,
            function: f,
            kind: kind,
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

    pub fn compute(mut self) -> IntegrationResult {
        let mut result = f64::NAN;
        let mut abserr = f64::NAN;

        let mut f = gsl_sys::gsl_function_struct {
            function: Some(rust_integrand::<F>),
            params: &mut self.function as *mut _ as *mut raw::c_void,
        };

        let err = match self.kind {
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
            }
        };

        IntegrationResult { value: result, abserr: abserr }
    }
}


impl IntegrationWorkspace {
    pub fn qagiu<'a, F>(&'a mut self, f: F, lower_bound: f64) -> IntegrationBuilder<'a, F> where F: FnMut(f64) -> f64 {
        IntegrationBuilder::new(self, f, Integrator::QAGIU, lower_bound, f64::NAN)
    }
}
