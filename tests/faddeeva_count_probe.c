/* Test-only linker instrumentation: no production counter or callback cost. */
#include <complex.h>
static int calls;
extern double complex __real_plume_faddeeva_z_value(double, double complex);
double complex __wrap_plume_faddeeva_z_value(double k, double complex z)
{
    ++calls;
    return __real_plume_faddeeva_z_value(k,z);
}
int plume_test_z_calls(void) { return calls; }
