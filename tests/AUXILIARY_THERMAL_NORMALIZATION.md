# Stage 3B.1: auxiliary thermal normalization

## Audit and mathematical decision (before prototype implementation)

The auxiliary approach is mathematically valid for the collisionless kinetic
susceptibility. The physical reference mass, charge, density, gyrofrequency,
Alfvén speed and inertial length remain fixed. No auxiliary population is needed.
The equations in `disp0` can be reused unchanged once the normalization-only
anisotropy `alphp` is separated from the physical anisotropy of species 1.

The audit covers `normalization.md`, the Stage 0/1, 2 and 3A reports,
`vars::sync_solver_state`, input reading and collision initialization in
`functions`, both complete susceptibility bodies, and the normalization/mixed
reference tests. `disp0` and `disp1` contain **no direct use of species 1's
temperature** outside its ordinary species response. In particular, neither
requires the internal `tau_s` for species 1 to equal unity. Their special use
of species 1 is its anisotropy in `alphp` and the cache predicate.

The full `disp0` audit identifies the following appearances of thermal state:

- `betap`, `vtp` occur in the species prefactors, drift normalization, final
  longitudinal drift-current correction and electromagnetic wave terms.
- `tau_s`, `alphp` enter the Bessel and plasma-dispersion arguments and six
  component prefactors. `disp_alph=spec(is)%alph_s` is a separate **physical**
  species anisotropy throughout the harmonic sums.
- `alphp=spec(1)%alph_s` and the Bessel-cache reference-anisotropy comparison
  are implementation couplings. They do not require the physical temperature
  of species 1 to define the normalization. Overwriting species 1's physical
  anisotropy would, however, change its response and is not permissible.
- Species mass, signed charge, density and Alfvén-normalized drift retain their
  physical meanings. Frequencies and gyrofrequencies need no conversion.
- `solver_*` aliasing is the global-state interface used by the legacy kernel;
  physical input state remains distinct. `sync_solver_state` is unsuitable for
  the cold-reference prototype and must not be called by it.
- `disp1` uses the same normalization combinations, with additional fixed
  `nu_ns/Omega_R` terms. Its algebra is likewise invariant if this *physical*
  collision rate is held fixed. Converting a Knudsen number as though the
  auxiliary length were a physical input length would change the problem.
  No collisional prototype or production collision change is made here.
- Diagnostics, FPC, scans and root workflows have additional thermal assumptions
  and are outside this isolated susceptibility experiment.

For a general positive normalization anisotropy `h`, take

\[
 B=b_*,\quad V=a\sqrt B,\quad K_j=\kappa_j\sqrt{Bh},\quad
 t_s=B/\theta_s,\quad \theta_s=\tau_s^{(d)},\quad a=v_A/c.
\]

`B,V,K,t` stand for the legacy solver's beta, thermal speed over c,
wavevectors and inverse-temperature ratios. The proposed isotropic auxiliary
scale is exactly `h=1`; physical `alpha_s` is never changed. The original
thermal normalization is recovered with `B=beta_R` and `h=alpha_R`.

Define the following invariant physical quantities:

\[
 g_s=\mu_s/Q_s,\quad A_s=D_s\mu_s/(Q_s^2a^2),\quad
 d_s=w-\kappa_\parallel u_s,\quad
 C_s=1/[\kappa_\parallel\sqrt{\mu_s\theta_s}],
\]
\[
 \lambda_s=\frac{\kappa_\perp^2 Q_s^2\alpha_s\theta_s}{2\mu_s},\quad
 \zeta_{ns}=C_s(d_s-ng_s).
\]

The cancellations are explicit:

\[
 B/V^2=1/a^2,\quad K_\parallel u_s/\sqrt{Bh}=\kappa_\parallel u_s,
\]
\[
 \frac{\sqrt{h t_s/\mu_s}}{K_\parallel}=C_s,\quad
 \frac{K_\perp^2 Q_s^2\alpha_s}{2\mu_s t_s h}=\lambda_s,
\]
\[
 \frac{v_{t\parallel,s}}c=V\sqrt{\mu_s/t_s}
       =a\sqrt{\mu_s\theta_s},\qquad
 v_{t\perp,s}=\sqrt{\alpha_s}\,v_{t\parallel,s}.
\]

Thus both thermal speeds, all resonances and Bessel arguments are physical.
The gyrofrequency sign is preserved. The causal `zet_in` branch sees the same
sign of parallel wavevector because `sqrt(Bh)>0`.

All six component prefactors in `disp0` reduce to

\[
 N_{xx}=N_{yy}=\frac{A_s}{w^2\lambda_s},\quad
 N_{xy}=-\frac{i A_s C_s}{w^2},\quad
 N_{xz}=\frac{D_s\kappa_\perp}{a^2w^2Q_s\lambda_s\kappa_\parallel},
\]
\[
 N_{yz}=\frac{iD_s\kappa_\perp}{a^2w^2Q_s\kappa_\parallel},\quad
 N_{zz}=\frac{2D_s}{a^2w^2Q_s^2\kappa_\parallel^2\alpha_s\theta_s}.
\]

The final longitudinal drift-current term reduces to

\[
 \Delta\chi_{zz}=\frac{2u_sD_s}{a^2w\alpha_sQ_s^2\kappa_\parallel\theta_s}.
\]

Every remaining term of the six harmonic sums is a polynomial in `w`, `d_s`,
`n*g_s`, `alpha_s`, `C_s`, `lambda_s`, the scaled Bessel functions and
`Z(zeta_ns)`. Consequently `chi_xx=N_xx H_xx`, `chi_yy=N_yy H_yy`,
`chi_xy=N_xy H_xy`, `chi_xz=N_xz H_xz`, `chi_yz=N_yz H_yz`, and
`chi_zz=N_zz H_zz+Delta chi_zz` are individually invariant, including the
zero harmonic. The remaining three elements obey the unchanged transpose
relations `chi_yx=-chi_xy`, `chi_zx=chi_xz`, `chi_zy=-chi_yz`, without conjugation.
This accounts for every complex element, not just the determinant.

For an explicit compact form of all elements, suppress the species subscript,
let `v=sqrt(mu*theta)`, and define

\[
 \Gamma_n=e^{-\lambda}I_n(\lambda),\quad
 H_n=\Gamma_n-\tfrac12(\Gamma_{n-1}+\Gamma_{n+1}),\quad Z_n=Z(\zeta_n),
\]
\[
 F_n=\alpha-1+[\alpha d+(1-\alpha)ng]C Z_n,
\]
\[
 B_n=\frac{\alpha(w-ng)-\kappa_\parallel u+ng
       +(w-ng)[\alpha d+(1-\alpha)ng]C Z_n}{\kappa_\parallel}.
\]

Here `Z_n` uses the causal continuation for the sign of `kpar`. With
`chi=(A/w^2)Y`, the equivalent harmonic-sum form used by the independent
NHDS/Stix reference is

\[
 Y_{xx}=\sum_n\frac{n^2\Gamma_n}{\lambda}F_n,\qquad
 Y_{xy}=-i\sum_n nH_nF_n,
\]
\[
 Y_{xz}=\frac{\kappa_\perp}{g\lambda}\sum_n n\Gamma_n B_n,\qquad
 Y_{yy}=\sum_n\left(\frac{n^2\Gamma_n}{\lambda}+2\lambda H_n\right)F_n,
\]
\[
 Y_{yz}=\frac{i\kappa_\perp}{g}\sum_n H_nB_n,\qquad
 Y_{zz}=\frac{2}{\kappa_\parallel v^2\alpha}
       \left[\sum_n(w-ng)\Gamma_n B_n+wu\right].
\]

The other elements are `Y_yx=-Y_xy`, `Y_zx=Y_xz`, `Y_zy=-Y_yz`.
All symbols are physical invariants; neither `B=b_*` nor `h=alpha_*` remains.
These equations document the cancellation and the independent oracle; they
are **not** a replacement kinetic implementation in Fortran. The existing
expanded paired-harmonic expressions and harmonic cutoff remain unchanged.

Even the legacy wave terms reduce to
`K_i*K_j/(w*V)^2/h = kappa_i*kappa_j/(a*w)^2`.
The prototype inspects individual susceptibilities only; it does not enable
new physical dispersion configurations or roots in production.

## Isolated prototype architecture

`src/auxiliary_normalization.f90` is an independently callable module, omitted
from the production Makefile's object list. Its interface separates native
inputs, auxiliary scales and derived species parameters:

```fortran
call make_auxiliary_scale(kpar_d, kperp_d, va_ref_c, b_star, scale, ierr)
call scale%export(beta_solver, vtp_solver, kpar_solver, kperp_solver, alpha_star)
call auxiliary_species(physical_species, scale, kinetic_species, ierr)
```

The private fields of `type(auxiliary_scale)` are populated by the constructor;
`alpha_star=1` is explicit on export. Neither operation accesses any global
plasma state or edits its inputs. Conversion copies **one** `specie` value,
changes only its inverse-temperature field and retains physical anisotropy,
mass, signed charge, density, drift and collision-rate field. The collisionless
experiment does not evaluate that collision rate. No full species-array copy
is needed in the adapter and no global auxiliary-temperature parameter exists.
The driver supplies arrays because that is the current kernel interface.

Adapter status codes are 0 success, 1 invalid parameters, 2 not a warm species
(exactly zero temperature), and 3 unrepresentable converted scale/ratio.
Negative/nonfinite temperatures are invalid. A cold species returns before
any temperature division or anisotropy access. Check status before exporting
or using derived data. Finite positive `b_*` and Alfvén speed are required;
used warm-species parameters must also be finite and physically valid. The
adapter checks representability of its own outputs, not of every later
legacy intermediate product. No numerical floors are introduced.

**No existing production source is modified.** The initial optional-argument
experiment in the large `disprels` compilation unit caused roundoff changes
in unrelated density/velocity diagnostics under GNU optimization, violating
the exact legacy regression. That edit was reverted in full. This is a
compiler/regression isolation constraint, not a failure of the mathematics.

`tests/auxiliary_kernel.py` instead extracts the current `disp0` body at test
build time. It changes only the routine name, adds a required normalization
anisotropy argument, binds `alphp` to that argument, and removes the legacy
Bessel-cache lookup. Bessel reuse stays false. It imports the existing `bessel`
and `zet_in` implementations and the existing harmonic cutoff. It does not
modify any susceptibility expression or keep a second checked-in copy of
the kinetic equations. Assertions and exact source anchors fail if the expected
structure changes. The generated module is compiled only into the test probe,
not the production executable, and has its own work/cache arrays.

`tests/auxiliary_probe.f90` manually constructs the physical test state and
bridges the adapter outputs into the existing `solver_*` interface. It does not
invoke `sync_solver_state` for auxiliary calls. It initializes solver scales
to NaN before each auxiliary conversion, keeps physical beta zero for a cold
reference, and tests cold anisotropy placeholders including NaN. Bitwise
snapshots verify that all physical globals and species fields remain unchanged.
Actual cold populations are skipped by the kinetic loop; their work-array rows
are zeroed only for the kernel's unused accumulation. No cold response, total
wave matrix, determinant or root is reported as a physical mixed calculation.
Only warm-species susceptibility elements are emitted.

For warm-reference comparison calls the probe can also invoke the **unmodified**
production `disp`, or the extracted body with the exact original thermal state
and original reference anisotropy. Returning to production after an auxiliary
call reproduces the earlier production susceptibility bit for bit, testing
cache separation. A process uses one fixed species count, matching the legacy
cache allocation assumptions. The adapter is stateless; the experimental
kernel binding is global-state based and is not a reentrant library interface.

## Automated validation

`test_auxiliary_normalization.py` builds production-flag and checked versions.
Both use compiler-promoted double precision (`nrtype::SP`, 64-bit real), not
`DP`, which becomes wider under the existing GNU flags. Adapter IEEE validation
is compiled with `-fno-fast-math`. The checked build uses
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fdefault-real-8`.
The production experiment uses `-O3 -funroll-loops -ffast-math -fdefault-real-8`
and links the normal production special-function objects.

Nine new tests cover:

- Fixed physical plasmas at `b_*=0.01, 0.1, 1, 10`: isotropic and anisotropic
  proton/electron systems; electron drift; alpha drift; three and four warm
  species; oblique wavevectors of both parallel signs; growing/damped complex
  evaluation frequencies; `v_A/c=0.001` and `0.05`.
- Auxiliary results against unmodified production kinetic susceptibilities
  for every warm species and every complex element. The exact original
  normalization (`B=beta_R`, `h=alpha_R`) and the isotropic auxiliary scale
  (`B=beta_R`, `h=1`) are both compared for an anisotropic reference.
- Cold protons with warm electrons, drifting electrons, anisotropic electrons,
  warm electrons plus warm alphas, and cold protons/cold electrons with warm
  drifting alphas. Every warm component is compared against the existing
  independent SciPy NHDS/Stix harmonic-sum implementation in `mixed_reference.py`.
  The first warm species' temperature is also tested as `b_*`.
- Physical species temperatures spanning `0.001` to `100` in one plasma,
  while remaining kinetic regardless of the auxiliary scale.
- Cold-reference anisotropy placeholders NaN, 0, -3, `1e100`, with identical
  results under floating-point traps; poisoned solver scales; bitwise unchanged
  physical state and rejection of invalid auxiliary parameters.
- Broad and extreme auxiliary-scale conditioning, failed arithmetic followed
  by a valid recomputation in the nontrapping build, and explicit traps in the
  checked build when the legacy arithmetic exceeds its representable range.
- The actual production input path still rejecting cold-reference/warm-species
  combinations, and `nm` confirming neither experimental module is linked into
  `plume.e`.

The warm production and scale comparisons use componentwise
`atol=1e-10, rtol=5e-9`. Independent NHDS/Stix comparisons use
`atol=1e-5, rtol=2e-7`, as appropriate to the large ordinary susceptibilities
and existing kinetic special-function approximations at `v_A/c=0.001`.
These are new comparisons; existing regression tolerances and expected values
are untouched. All nine elements are compared; no determinant-only proxy is used.

## Numerical conditioning and measured results

The measured discrepancy is `abs(actual-reference)/max(1,abs(reference))`,
maximized over every complex element in both builds:

| Comparison | Maximum discrepancy |
| --- | --- |
| Standard auxiliary scales for a fixed physical plasma | `4.073e-9` |
| Auxiliary versus unmodified production kinetic calculation | `2.036e-9` |
| Cold-reference prototype versus independent NHDS/Stix reference | `6.223e-9` |
| Wide-scale stress comparison versus `b_*=1` | `2.548e-9` |

The scale dependence is mathematically absent but not bitwise absent from
floating-point arithmetic. For example, in the cold-proton/warm-electron/warm-
alpha stress case, changing `b_*` from 1 to `1e-12` changes the electron `chi_xx`
by `1.145e-5` in absolute units, or `2.548e-9` relative to its magnitude of
about `4494`. Its physical Bessel argument is only `1.089e-7`. This is consistent
with cancellation in the unchanged anisotropic electron paired-harmonic sums:
order-unity anisotropy terms cancel to produce a small net transverse response.
The alpha response in the same comparison differs by at most `2.884e-14`
scaled. Different auxiliary scales introduce rounding into the same physical
Bessel/resonance arguments and prefactors; no physical temperature has changed.

The stress state uses cold protons, electrons with
`theta_e=0.2`, `alpha_e=0.8`, `mu_e=1836`, `u_e=0.1`, and alphas with
`theta_alpha=0.03`, `alpha_alpha=1.2`, `mu_alpha=0.25`, `Q_alpha=0.5`,
`D_alpha=0.05`, `u_alpha=0.3`. Evaluation is at
`w=0.13+0.03i`, `(kpar*d_R,kperp*d_R)=(0.06,0.05)`, `v_A/c=0.001`.

| Auxiliary scales / build | Observation |
| --- | --- |
| `1e-70, 1e-40, 1e-12, 1e-6, 0.01, 0.1, 1, 10, 1e6, 1e12, 1e40, 1e70`, both builds | Finite tensors; pass the scale-comparison tolerance. |
| `1e-80`, `1e80`, optimized production flags | Nonfinite intermediate arithmetic; probe reports status 4, then a valid-scale call recovers. |
| `1e-80`, `1e80`, checked build | Finite; discrepancies `3.577e-14`, `2.104e-13`. |
| `1e-160`, `1e160`, optimized flags | Nonfinite intermediate arithmetic; subsequent valid-scale evaluation recovers. |
| `1e-160`, `1e160`, checked build | Explicit floating-point traps in the legacy arithmetic. |

These boundaries are measured for this state/compiler, not universal allowable
ranges. Even when the adapter outputs are finite, unsimplified kinetic products
such as `B*t_s` and `V^2*Kpar^2`, and complex-division intermediates, can overflow
or underflow before their mathematical cancellation. Optimization and fast-math
change that range; the checked build's range can be wider for some expressions.
The tests allow a future compiler to evaluate an extreme scale successfully,
but then require its tensor to agree, rather than pinning a failure threshold.
No production arithmetic, damping, temperature floor or auxiliary-scale
selection optimization is introduced to hide these limitations.

`lambda_s` and `zeta_ns` themselves are physical, independent of `b_*`.
Thus changing `b_*` cannot cure a genuinely large Bessel argument, inadequate
harmonic truncation, a large/resonant plasma-dispersion argument, or cancellation
in a physically small susceptibility. The test uses PLUME's existing 128-harmonic
cutoff and special functions. The prototype inherits their conditioning limits.

A first warm species' temperature is a reasonable **provisional** future scale
for ordinary representable parameter ranges: it gives that species `t_s=1`,
avoids tying the scale to the cold reference, and passed the cold-reference
comparisons. The five-decade temperature test also passes. This is not a
universal conditioning guarantee for arbitrary temperature spreads. Selection
policy and broader extreme-parameter handling remain future work.

## Regression results and Stage 3B.2 recommendations

With GNU Fortran 11.4.0, `make -j2` succeeds and the requested command

```sh
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

passes **59 tests**: all 50 prior tests and nine auxiliary-normalization tests.
The pinned legacy comparison is **exact**, with no changed expectations.
Existing mixed and fully cold benchmarks pass unchanged. The ordinary
normalization suite retains its `4.668e-13` maximum equivalence discrepancy.
Python compilation and `git diff --check` pass. Intel Fortran and rendered
FORD documentation were not run; those tools are unavailable locally.

Only five files are added by Stage 3B.1:

| File | Purpose |
| --- | --- |
| `src/auxiliary_normalization.f90` | Independent scalar/one-species conversion API; not in production build. |
| `tests/auxiliary_kernel.py` | Extract and bind the existing kinetic body for isolated tests, without maintaining duplicate equations. |
| `tests/auxiliary_probe.f90` | Direct susceptibility probe, physical-state preservation, production comparison and status reporting. |
| `tests/test_auxiliary_normalization.py` | Nine discoverable validation/conditioning/isolation tests. |
| `tests/AUXILIARY_THERMAL_NORMALIZATION.md` | Audit, transformations, implementation and measured results. |

No existing production file, Makefile, input definition, root method, scan,
cold kernel, diagnostic, heating or FPC implementation is changed. Unrelated
workspace edits and pre-existing untracked files are preserved.

For Stage 3B.2:

1. Retain all physical reference definitions and classify species before the
   thermal adapter. Reuse the explicit auxiliary conversion API only for warm
   components of a cold-reference state; evaluate cold responses in native
   inertial units and keep the existing supported paths unchanged.
2. Separate normalization anisotropy from species anisotropy in a production
   kinetic-core interface. Isolate compilation carefully and require exact
   legacy regression: mathematical equivalence alone did not preserve optimizer
   decisions in the large existing module. The test-generated kernel is an
   experimental access mechanism, not the proposed production dispatch.
3. Make cache keys reflect the actual Bessel argument or explicit normalization
   state, and handle species-count allocation explicitly. Do not substitute a
   fictitious particle species or overwrite physical `alpha_R`.
4. Sum the independently evaluated cold and kinetic susceptibilities before
   one native electromagnetic/vacuum assembly. Continue to use the existing
   root finder and status mechanisms after separate integration validation.
5. Exact parallel kinetic propagation can use the existing native thermal-speed
   limit, independent of reference temperature, but needs dedicated integration
   tests. This general-`disp0` prototype deliberately requires nonzero `kpar`
   and `kperp`; it does not implement their singular kinetic limits. Warm exact
   perpendicular propagation, collisions, scans, eigenfunctions, heating and
   FPC remain outside this task.

The auxiliary approach is valid and the kinetic equations are reusable.
Remaining obstacles are explicit normalization/cache interfaces, floating-point
conditioning and production regression isolation, not a need for a native
inertial-variable rewrite. Stage 3B.2 is not activated here.
