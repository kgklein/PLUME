/* PLUME ISO_C_BINDING bridge. Numerical algorithm: unmodified MIT Faddeeva.
 * Keep the whole normalized/signed evaluation in an IEEE environment: GNU's
 * legacy -ffast-math link startup enables both FTZ and DAZ. Fortran's gradual
 * underflow control alone does not disable DAZ on this toolchain. */
#include <complex.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "Faddeeva.h"

static double complex pair(double re, double im)
{
    /* C99 defines double complex with the representation of double[2], real
       part first. Copying also constructs NaNs without the 0*Inf pitfall. */
    const double parts[2] = {re, im};
    double complex z;
    memcpy(&z, parts, sizeof z);
    return z;
}

double complex plume_faddeeva_z(double kpar, double complex z, int *status)
{
    const double sqrt_pi = 1.77245385090551602729816748334114518;
    double complex result = pair(NAN, NAN), w;
    fenv_t saved;
    int rounding;
    *status = 4; /* numerical failure */
    if (fegetenv(&saved) != 0) return result;
    rounding = fegetround();
    if (fesetenv(FE_DFL_ENV) != 0) return result;
    if (rounding != -1 && fesetround(rounding) != 0) goto finish;
    *status = 1; /* invalid input */
    if (!isfinite(kpar) || !isfinite(creal(z)) || !isfinite(cimag(z))) goto finish;
    if (kpar == 0.0) { *status = 2; goto finish; }
    if (kpar < 0.0) z = -z;
    w = Faddeeva_w(z, 0.0); /* upstream's full machine accuracy */
    if (isnan(creal(w)) || isnan(cimag(w))) { *status = 4; goto finish; }
    if (!isfinite(creal(w)) || !isfinite(cimag(w)) ||
        fabs(creal(w)) > DBL_MAX/sqrt_pi || fabs(cimag(w)) > DBL_MAX/sqrt_pi) {
        *status = 3; /* output exceeds binary64 component range */
        goto finish;
    }
    result = pair(-sqrt_pi*cimag(w), sqrt_pi*creal(w));
    if (kpar < 0.0) result = -result;
    if (!isfinite(creal(result)) || !isfinite(cimag(result))) {
        result = pair(NAN, NAN); *status = 4;
    } else *status = 0;
finish:
    /* fesetenv restores the caller's modes AND prior flags, not this call's
       intermediate exceptions. It does not raise them as feupdateenv would. */
    if (fesetenv(&saved) != 0) { *status = 4; result = pair(NAN, NAN); }
    return result;
}

/* Side-effect-free value interface for the pure Fortran kinetic callback.
 * Failure is encoded as a quiet NaN and classified by the evaluation guard. */
double complex plume_faddeeva_z_value(double kpar, double complex z)
{
    int status;
    return plume_faddeeva_z(kpar,z,&status);
}
