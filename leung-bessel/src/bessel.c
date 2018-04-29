/* A fast Bessel function calculator, useful for evaluating Bessel functions
 * J(n, z) and J'(n, z) for large n.
 *
 * The original implementation is by Po Kin Leung and is described in Leung,
 * Gammie, and Noble (2011, DOI:10.1088/0004-637X/737/1/21) and was made
 * public in their "harmony" code. This version is the deriviative from its
 * descendant, "symphony", described by Pandya, Zhang, Chandra, and Gammie
 * (2016, DOI:10.3847/0004-637X/822/1/34). Specifically, the version in Git
 * commit 41a42e545b003e7df71ad286b1d4bba2fa329702 of
 * https://github.com/AFD-Illinois/symphony/. It has been tidied up by Peter
 * K. G. Williams.
 *
 * This version is distributed under the GPL version 3 and is Copyright 2011 Po
 * Kin Leung, 2016 Alex Pandya, and 2017 Peter Williams.
 */

#include <float.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

/* Compute `f_factor * exp(f_exp)` with some extra precision and safety. */
static inline double
exp_factor(const double f_factor, const double f_exp)
{
    double fabs_exp;

    if (f_factor == 0.)
        return 0.;

    fabs_exp = fabs(f_exp);

    /* If small, use a 9th-order Taylor expansion of the exponent. */
    if (fabs_exp < 1e-3) {
        const double x = f_exp;
        return f_factor *
            (1 + ((40320 + (20160 + (6720 + (1680 + (336 + (56 + (8 + x) * x) * x) * x) * x) * x) * x) * x / 40320.));
    }

    /* We want to prevent overflows, calculate exponential wisely if we are close: */
    if (fabs_exp > 690.) {
        const double sign_f = (f_factor < 0) ? -1. : 1.;
        const double log_f = log(fabs(f_factor));

        if (log_f * f_exp < 0.)
            return sign_f * exp(log_f + f_exp);
        return f_factor * exp(f_exp);
    }

    /* Standard expression: */
    return f_factor * exp(f_exp);
}


/* Meissel's "second" expansion, specified to higher order by Chishtie et al.
 * 2005. Good for x >> n.
 */
static inline double
BesselJ_Meissel_Second(const double n, const double x)
{
    const double z = x / n;
    const double eps = (x - n) / n;
    const double Z = sqrt(eps * (1 + z));
    const double U = 1. / (n * Z * Z * Z);
    const double t1 = z * z;
    const double t2 = U * U;

    // P_n sum (exponent = -Psum, so this is -Psum)
    const double exp_val =
        (t1 * t2 * (-3072 - 768 * t1 + (3072 + (27648 + (22272 + 1248 * t1) * t1) * t1 +
        (-3072 + (-165120 + (-952576 + (-1119552 + (-271488 - 6592 * t1) * t1) * t1) * t1) *
         t1 + (3072 + (744960 + (15287808 + (72179904 + (102842688 + (45756144 + (5297808 +
        71391 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t2) * t2) * t2)) / 0.12288e5;

    // Part 1 of Q_n sum
    const double Qt = n * (Z - acosl(n / x));

    // Part 2 of Q_n sum
    const double Qsum =
        -(U * (860160 + 1290240 * t1 + (28672 + (-2709504 + (-6547968 - 672000 * t1) * t1) *
        t1 + (8192 + (2519040 + (60518400 + (151828480 + (61254720 + 2163168 * t1) * t1) * t1) *
        t1) * t1 + (-6144 + (2644992 + (299351808 + (3405435264 + (8653594320 + (5897669400 +
        (954875250 + 16907985 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t2) * t2) * t2)) /
        0.10321920e8;

    const double factor = sqrt(2 / (M_PI * n * Z)) * cosl(Qsum + Qt - M_PI_4);

    return exp_factor(factor, exp_val);
}


/* Meissel's "first" expansion, specified to higher order by Chishtie et al.
 * 2005. Good for x << n.
 */
static inline double
BesselJ_Meissel_First(const double n, const double x)
{
    double exp_val;

    const double z = x / n;
    const double eps = (n - x) / n;
    const double Z = sqrt(eps * (1 + z));
    const double ninv = 1. / n;
    const double U = 1. / (n * Z * Z * Z);
    const double t1 = z * z;
    const double t2 = ninv * ninv;

    // Part 1 of V_n sum
    const double Vsum1 =
        (U * (860160 + 1290240 * t1 + ((-2580480 - 645120 * t1) * t1 + (-28672 + (2709504 +
        (6547968 + 672000 * t1) * t1) * t1 + ((-2580480 + (-23224320 + (-18708480 -
        1048320 * t1) * t1) * t1) * t1 + (-8192 + (-2519040 + (-60518400 + (-151828480 +
        (-61254720 - 2163168 * t1) * t1) * t1) * t1) * t1 + ((2580480 + (138700800 +
        (800163840 + (940423680 + (228049920 + 5537280 * t1) * t1) * t1) * t1) * t1) * t1 +
        (6144 + (-2644992 + (-299351808 + (-3405435264 + (-8653594320 + (-5897669400 +
        (-954875250 - 16907985 * t1) * t1) * t1) * t1) * t1) * t1) * t1 + (2580480 +
        (625766400 + (12841758720 + (60631119360 + (86387857920 + (38435160960 +
        (4450158720 + 59968440 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t1 * U) * U) *
        U) * U) * U) * U) * U)) / 0.10321920e8;

    // Part 2 of V_n sum
    const double Vsum2 = -(ninv * (420 + (-14 + (-4 + 3 * t2) * t2) * t2)) / 0.5040e4;

    // I substitute Gamma(n+1) with (n+1)*Gamma(n) in the denominator:
    const double factor = 1. / ((n + 1.) * sqrt(Z));

    if (eps < 1e-4 && n > 1e3) {
        const double t3 = t2 * t2;
        // O(1/n^k) part of -log(gamma(n)):
        const double loggamma_exp = (ninv * (-420 + 14 * t2 - 4 * t3 + 3 * t3 * t2)) / 0.5040e4;
        const double exp2 =
            -n * sqrt(2.*eps) * eps * (0.984023040e9 + (0.442810368e9 + (0.303114240e9 +
            (0.233192960e9 + (0.190139040e9 + (0.160692840e9 + 0.139204065e9 * eps) * eps) *
            eps) * eps) * eps) * eps) / 0.1476034560e10;

        exp_val = 0.5 * log(0.5 * n / M_PI) + loggamma_exp + exp2 - Vsum1 - Vsum2;
    } else {
        double invZp1;

        // 1 / (1 + Z) to eigth order:
        if (Z < 1.e-3)
            invZp1 = 1 + (-1 + (1 + (-1 + (1 + (-1 + (1 - Z) * Z) * Z) * Z) * Z) * Z) * Z;
        else
            invZp1 = 1. / (1. + Z);

        exp_val = n * (log(x * invZp1) - (1 - Z)) - Vsum1 - Vsum2 - lgamma(n);
    }

    return exp_factor(factor, exp_val);
}


/* The Debye "epsilon" expansion. Good for in the limit x -> n from either
 * side but only valid very close to the limit; expansion is in the variable
 * "eps = x - n". Because of the order of the epsilon used, this routine is
 * limited to x < 1e55.
 */
#define BESSEL_EPSILON_ORDER 16

static inline double
BesselJ_Debye_Eps_Exp(const double n, const double x)
{
    static int first_time = 1;
    static double At[BESSEL_EPSILON_ORDER];

    if (x > 1.e55)
        return NAN;

    if (first_time) {
        int m;

        for (m = 0; m < BESSEL_EPSILON_ORDER; m++) {
            double M = (m + 1) / 3.;
            At[m] = sin(M_PI * M) * pow(6., M) * exp(lgamma(M));
        }

        first_time = 0;
    }

    const double ez = x - n;
    const double z = pow(x, 1. / 3.);
    const double t3 = z * z;
    const double t4 = x * z;
    const double t10 = t4 * t4;
    const double t38 = 810485676000000 * At[4] * t3;
    const double t44 = x;
    const double t70 = At[7] * t3;
    const double t93 = At[10] * t3;
    const double t107 = At[13] * t3;
    const double t114 = ez * ez;
    const double t146 = t10 * t10;
    const double t149 =
        ((-5403237840000 * At[6] + (69470200800000 * At[4] + 19451656224000000 *
        At[0] * t4) * t3) * t10 * z + (-401283384 * At[15] + (4707059994 * At[13] +
        (3012121710 * At[12] + (-36011689560 * At[10] + (8027667648000 * At[9] +
        (-8027667648000 * At[7] + (-1296777081600000 * At[3] + 19451656224000000 *
        At[1] * t3) * t4) * t3) * z) * t3) * z) * t3 + ((-41423013450 * At[12] +
        (484040056500 * At[10] + (67540473000000 * At[6] - t38) * t4) * t3) * t44 +
        (1748257220 * At[15] + (-19964735910 * At[13] + (-2594411820000 * At[9] +
        (29331862560000 * At[7] + 3241942704000000 * At[3] * t4) * t3) * t4) * t3 +
        ((78248884350 * At[12] + (-860873013000 * At[10] + (-94556662200000 * At[6] +
        t38) * t4) * t3) * t44 + (-1938419560 * At[15] + (20997160275 * At[13] +
        (2283511230000 * At[9] - 21612951360000 * t70) * t4) * t3 + ((-47153256150 *
        At[12] + (459918459000 * At[10] + 27016189200000 * At[6] * t4) * t3) * t44 +
        (849093050 * At[15] + (-8397889500 * At[13] + (-643242600000 * At[9] +
        3859455600000 * t70) * t4) * t3 + ((11448186750 * At[12] - 88445857500 * t93) *
        t44 + (-173573400 * At[15] + (1474097625 * At[13] + 53603550000 * At[9] * t4) *
        t3 + ((-1161410250 * At[12] + 5360355000 * t93) * t44 + (-113704500 * t107 +
        17481100 * At[15] + (40608750 * At[12] * t44 + (-833000 * At[15] + 3123750 *
        t107 + 14875 * At[15] * t114) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) *
        ez) * ez) * ez) * ez) * ez) / (M_PI * t146 * 0.58354968672000000e17);

    return t149;
}


/* The swiss-army-knife version. Requires x and n >= 0. Only especially good
 * for n > 100. The approximations are based on those given in Chishtie et a.
 * 2005.
 *
 * The parameters below are lines that deliminate between various
 * approximations to J_n(x), found empirically for n = 100..1e7.
 *
 * For our particular application, the approximations are only used for n >
 * 30; below this value, we use GSL, which only accepts integer `n` values.
 * Therefore non-integer `n` values below 30 will give NAN back.
 */

#define SLOPE2 -6.656260931106707801e-01
#define SLOPE3 -6.543033585865805080e-01
#define INTERCEPT2 2.563324856985127465e-01
#define INTERCEPT3 2.720161927383055733e-01
#define N_JN 30.

double
my_Bessel_J(const double n, const double x)
{
    double logn;

    if (!(n >= 0 && x >= 0))
        /* The above definition catches NAN inputs */
        return NAN;

    if (n < N_JN) {
        int n_int = (int) n;

        if (n_int != n)
            return NAN;

        return gsl_sf_bessel_Jn(n_int, x);
    }

    logn = log10(n);

    if (x < n) {
        double y = log10((n - x) / n);

        if (x != n && y > SLOPE2 * logn + INTERCEPT2)
            return BesselJ_Meissel_First(n, x);

        return BesselJ_Debye_Eps_Exp(n, x);
    } else {
        if (x != n) {
            double y = log10((x - n) / x);

            if (y > SLOPE3 * logn + INTERCEPT3)
                return BesselJ_Meissel_Second(n, x);
        }

        return BesselJ_Debye_Eps_Exp(n, x);
    }
}


/* The swiss-army-knife Bessel derivative. It is based on the recurrence relation:
 *
 *   J_n'(x) = -J_{n+1}(x) + J_n(x)*(n/x)
 *
 * It would also be possible to use:
 *
 *   2 J_n'(x) = J_{n-1}(x) + J_{n+1}(x)
 *
 * But the current choice lets us avoid special-casing n = 0.
 */

double
my_Bessel_dJ(const double n, const double x)
{
    const double jn = my_Bessel_J(n, x);
    const double jnp1 = my_Bessel_J(n + 1, x);

    if(x == 0.) {
        /* J_n(0) = 0 for n >= 1; the recurrence relation gives a zero derivative for n >= 2 */
        if (n >= 2.)
            return 0.;

        /* d(J_0(z))/dz = -J_1(z) */
        if (n == 0.)
            return -jnp1;

        return n * jn / DBL_MIN - jnp1;
    }

    return n * jn / x - jnp1;
}


/* ~Improved version
 *
 * TODO: The second Meissel approximation seems to have real problems for n >~
 * 10^9. At least I see them where we try to overlap it with Debye.
 */

const double MINUS_ETA_A_INTERCEPT = 0.174857;
const double MINUS_ETA_B_INTERCEPT = 0.295966;
const double PLUS_ETA_A_INTERCEPT = 0.151550;
const double PLUS_ETA_B_INTERCEPT = 0.438914;

double
pkgw_bessel_j(const double n, const double x)
{
    double logn;

    if (!(n >= 0 && x >= 0))
        /* The above definition catches NAN inputs */
        return NAN;

    if (n < N_JN) {
        int n_int = (int) n;

        if (n_int != n)
            return NAN;

        return gsl_sf_bessel_Jn(n_int, x);
    }

    if (x == n)
        return BesselJ_Debye_Eps_Exp(n, x);

    logn = log10(n);

    if (x < n) {
        const double eta = log10((n - x) / n); /* called "y" by Leung */
        const double eta_thresh_lo = -0.6666666 * logn + MINUS_ETA_A_INTERCEPT;
        const double eta_thresh_hi = -0.6666666 * logn + MINUS_ETA_B_INTERCEPT;

        if (eta < eta_thresh_lo)
            return BesselJ_Debye_Eps_Exp(n, x);

        if (eta > eta_thresh_hi)
            return BesselJ_Meissel_First(n, x);

        {
            const double debye = BesselJ_Debye_Eps_Exp(n, x);
            const double meissel1 = BesselJ_Meissel_First(n, x);
            const double pos = (eta - eta_thresh_lo) / (MINUS_ETA_B_INTERCEPT - MINUS_ETA_A_INTERCEPT);
            return debye * (1 - pos) + meissel1 * pos;
        }
    } else {
        const double eta = log10((x - n) / x);
        const double eta_thresh_lo = -0.6666666 * logn + PLUS_ETA_A_INTERCEPT;
        const double eta_thresh_hi = -0.6666666 * logn + PLUS_ETA_B_INTERCEPT;

        if (eta < eta_thresh_lo)
            return BesselJ_Debye_Eps_Exp(n, x);

        if (eta > eta_thresh_hi)
            return BesselJ_Meissel_Second(n, x);

        {
            const double debye = BesselJ_Debye_Eps_Exp(n, x);
            const double meissel2 = BesselJ_Meissel_Second(n, x);
            const double pos = (eta - eta_thresh_lo) / (PLUS_ETA_B_INTERCEPT - PLUS_ETA_A_INTERCEPT);
            return debye * (1 - pos) + meissel2 * pos;
        }
    }
}

double
pkgw_bessel_dj(const double n, const double x)
{
    const double jn = pkgw_bessel_j(n, x);
    const double jnp1 = pkgw_bessel_j(n + 1, x);

    if(x == 0.) {
        if (n >= 2.)
            return 0.;

        if (n == 0.)
            return -jnp1;

        return n * jn / DBL_MIN - jnp1;
    }

    return n * jn / x - jnp1;
}
