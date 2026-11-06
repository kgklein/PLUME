<a id="electromagnetic-eigenfunctions-and-species-moments"></a>

# Electromagnetic eigenfunctions and species moments

[Manual home](./) · [Normalization](normalization.md) · [Output columns](output.md)

<a id="field-eigenfunction"></a>

## Field eigenfunction

For the species susceptibilities chi_s, PLUME forms epsilon=I+sum_s chi_s and
the complete electromagnetic matrix
\(\Lambda=\epsilon+\mathbf n\mathbf n^T-n^2I\), with
\(\mathbf n=\boldsymbol\kappa/(a w)\). The definitions of w, kappa and a are
on [Normalization](normalization.md). A mode has det(Lambda)=0 and
\(\Lambda\mathbf E=0\). This uses **all** cold and kinetic populations;
species-resolved tensors are retained, not replaced by a fitted total response.

The historical solver constructs component ratios using two-by-two minors
(cofactor elimination), with **Ex=1+0i**. It is not a general eigenvalue solver.
The cold-containing field module uses the same algebra, row scaling and equivalent
cofactor pivots to handle cases such as exact-parallel transverse polarization.
It tests the null residual and rejects unresolved rank/gauge cases. The legacy
all-kinetic path retains its arithmetic for exact comparisons. An independent
SVD is a validation reference, not the production field algorithm.

E=deltaE/deltaEx fixes the otherwise arbitrary overall complex amplitude and
phase. The relative complex component phases describe polarization. Compare
component ratios or vectors aligned by an overall phase; unaligned raw vectors
from another eigensolver need not match. Coordinate order is x,y,z with B0 along
positive z and k=(kperp,0,kpar). Interpret handedness using these axes and the
exp(i k.r-i omega t) convention, not an unstated observer convention.

<a id="magnetic-conversion-and-units"></a>

### Magnetic conversion and units

In Gaussian units physical Faraday's law is
\(\omega\,\delta\mathbf B=c\mathbf k\times\delta\mathbf E\).
For B=deltaB/deltaEx, the output relationship is

\[
 w\mathbf B=\frac{\boldsymbol\kappa\times\mathbf E}{a}
 =\frac{\mathbf K\times\mathbf E}
        {\mathrm{vtp}_{\rho}\sqrt{\alpha_R}}.
\]

The last form applies only to a finite thermal reference with rho_ref inputs;
K=k rho_R and vtp_rho=v_tparallel,R/c. The first form applies directly to d_ref,
including cold references. The auxiliary thermal scale cancels; it is not a new
physical E/B amplitude scale. B here is not normalized to B0. If observations
use SI, convert physical fields consistently before comparing electric/magnetic
amplitudes.

For E=(1,Ey,Ez) and k_y=0, Bx=-kappa_parallel Ey/(a w),
By=(kappa_parallel-kappa_perp Ez)/(a w), Bz=kappa_perp Ey/(a w).
At complex w these are complex field ratios. Heating energy uses a distinct
real-frequency magnetic reconstruction, discussed on [Heating](heating.md).

<a id="ex1-is-a-restriction"></a>

### Ex=1 is a restriction

A mode with Ex=0 cannot be represented by this gauge. Examples include an
exactly perpendicular cold ordinary mode and certain longitudinal/degenerate
polarizations. A dispersion root can be valid while its Ex=1 field is unavailable.
PLUME does not introduce another field normalization or a small Ex floor.
Off-mode matrices, unresolved nullvectors and nonfinite evaluations also yield
unavailable fields. Do not mistake NaN for a physical zero polarization.

For cold-containing fields, status meanings are: 0 success; 1 invalid
input/dispersion; 2 unresolved rank; 3 unavailable Ex gauge; 4 no acceptable
nullvector; 5 nonfinite arithmetic. The field module accepts a maximum
row-normalized null residual of 1e-8 and uses a 64-epsilon gauge/rank check.
These checks are not guarantees of root accuracy or independent branch identity.
A legacy kinetic output should still be checked for finite components and an
independent matrix residual where accuracy matters.

<a id="species-density-and-velocity-normalization"></a>

## Species density and velocity normalization

Let e0=deltaEx/B0, \(\mathbf v_s=\delta\mathbf U_s/(c e_0)\), and
\(\eta_s=(\delta n_s/n_{0s})/e_0\). These are exactly the complex velocity and
density output slots. They are not velocity divided by thermal speed or v_A,
and density is relative to **that species'** n0s, not n_R. The same normalization
applies to rho_ref and d_ref, warm and cold references, and either model.

<a id="kinetic-currentcontinuity-method"></a>

### Kinetic current/continuity method

The existing kinetic moment calculation obtains the normalized current from
that population's ordinary susceptibility at the accepted **complex mode** w:

\[
 \mathbf F_s=-i w a^2\frac{Q_s}{D_s}\chi_s\mathbf E,
 \qquad
 \eta_s=\frac{\boldsymbol\kappa\cdot\mathbf F_s}{a w},\qquad
 \mathbf v_s=\mathbf F_s-a u_s\eta_s\hat{\mathbf z}.
\]

Here Q_s is the signed input charge ratio, D_s the density ratio and u_s the
Alfvén-normalized parallel drift. F_s includes the equilibrium-density current.
This compact current removal/continuity form is used in the mixed path and is
algebraically equivalent to the historical kinetic formulas where their
intermediate divisions are regular. It reuses each kinetic susceptibility,
not a duplicate harmonic moment solver. The all-kinetic path preserves its
historical arithmetic. For D_s=0 relative density/bulk moments of an absent
population are undefined and marked unavailable on guarded paths.

<a id="cold-fluid-momentum-and-continuity"></a>

### Cold-fluid momentum and continuity

Define \(\delta_s=w-\kappa_\parallel u_s\) and \(g_s=\mu_s/Q_s\).
The independently normalized cold equations are

\[
 -i\delta_s\mathbf v_s=g_s\left[
 \mathbf E+\mathbf v_s\times\hat{\mathbf z}
          +a u_s\hat{\mathbf z}\times\mathbf B\right],\qquad
 \delta_s\eta_s=\frac{\boldsymbol\kappa\cdot\mathbf v_s}{a}.
\]

The magnetic drift force is essential for oblique waves. PLUME solves the
three-component momentum matrix and then continuity, without thermal speeds,
Z or Bessel sums. No cold anisotropy placeholder enters. The signs follow the
same Fourier convention and signed charge as the susceptibility kernel.
Exact delta=0 and delta=+/-g are poles; no frequency offset or damping hides them.
Nearby finite complex frequencies are evaluated directly.

The complete species current is

\[
 \delta\mathbf J_s=q_s[n_{0s}\delta\mathbf U_s
                 +\delta n_s U_s\hat{\mathbf z}],
 \qquad
 \frac{D_s}{Q_s}[\mathbf v_s+a u_s\eta_s\hat{\mathbf z}]
       =-i w a^2\chi_s\mathbf E.
\]

Dropping the density-drift term breaks consistency even if the transverse
velocities look plausible. The cold moment kernel is independently validated
against this susceptibility current, momentum and continuity. Individual
species indices remain those of the input; neither electrons nor cold species
are assumed to occupy particular slots.

<a id="requests-failure-isolation-and-examples"></a>

## Requests, failure isolation and examples

Serial scans use `eigen=.true.` in the relevant scan_input_N; two-dimensional
grids use scan_input_2 for eigen/heating flags. Direct option=-1 attempts
available diagnostics at the supplied frequency. option=1 with nscan=0 does
not write a standalone eigenfunction table. Map eigen/tensor requests are
ignored with warnings, not converted into diagnostic files.

Cold-containing dispatch calculates fields, then each species' warm or cold
moment. A failed species moment leaves that block NaN without overwriting valid
E/B or other species. A failed field leaves all associated moments unavailable.
Moment statuses: 0 valid; 1 input/classification; 2 zero w; 3 cold Doppler pole;
4 cold cyclotron pole; 5 nonfinite arithmetic; 6 absent population. Diagnostics
must not contaminate subsequent scan points. Corrected heating has separately
isolated state and cannot invalidate accepted roots/moments on failure.

```bash
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
./plume.e inputs/example/mixed_fields_cold_reference.in
./plume.e inputs/example/cold_fields_parallel.in
./plume.e inputs/example/mixed_moments_multispecies.in
```

These bounded examples write E/B plus kinetic/fluid moments in existing mode
columns; see [Examples](examples.md) for their actual inputs and benchmarks.
The cold parallel field example's first written root is approximately
0.05820986 in units of Omega_R, with imaginary part about -1.3e-19 (numerical
roundoff about a real cold branch). It is a transverse circular mode; no
kinetic damping should be inferred from that tiny imaginary value.

For an independent check, reconstruct Lambda from species chi and compare
\(\|\Lambda E\|/(\|\Lambda\|\|E\|)\), or the stricter row-scaled residual used by
the field module. Verify Faraday's law and align any independent nullvector to
Ex=1. A small raw determinant alone is insufficient when the matrix is badly
scaled. [Output](output.md) gives the column order; pressure perturbations,
new gauge choices and cold kinetic heating/FPC are not implemented.
