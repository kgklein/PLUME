<a id="developer-apis-verification-and-historical-reports"></a>

# Developer APIs, verification and historical reports

[Manual home](./) · [Numerics](numerics.md) · [Input contract](input.md)

The scientific user entry points are the manual and input workflows, not a list
of generic internal procedures. FORD-generated API documentation remains useful
for developers tracing exact signatures and source comments. It is not an
alternative authoritative definition when historical comments conflict with
current equations, output writers or tests.

<a id="dependency-map"></a>

## Dependency map

```
plume program -> functions: parse indexed namelists and validate
              -> map_search / refine_guess / scan loops
                 -> disp: composition/normalization dispatch
                    -> kinetic or auxiliary kinetic response
                    -> cold susceptibility
                    -> species susc + total lam + determinant
                 -> electromagnetic fields and species moments
                 -> corrected local heating/energy snapshots -> companion writer
                 -> mode/tensor writers
```

Species chi uses ordinary dimensionless susceptibility; lambda includes vacuum
and electromagnetic terms once. `solver_*` are derived thermal coordinates,
not authoritative physical inputs. Cold-reference and fully cold dispatch must
not construct a physical zero-temperature gyroradius. Heating derivatives must
not leave global state at the real or perturbed frequency rather than the
accepted complex mode.

<a id="api-entry-points"></a>

## API entry points

FORD links below resolve in a built current-source site. They retain the existing
modules/procedures rather than introducing a second documentation API.

| Interface | Purpose |
|---|---|
| [[disprels(module)]] / [[disp]] | Public production determinant callback, ordinary species susc and total lam |
| [[functions(module)]] / [[read_in_params]] | Input parser and indexed-group setup |
| [[vars(module)]] | Physical species/scanner types, input and derived state |
| [[calc_chi_cold]] | Explicit single-species cold response, native inertial coordinates, status codes |
| [[calc_warm_auxiliary]] | Warm response in explicit auxiliary coordinates, physical inputs read-only |
| [[calc_em_eigen]] / [[calc_em_fields]] | Existing Ex=1 field method with explicit physical normalization |
| [[calc_cold_moments]] / [[calc_kinetic_moments]] | Explicit species fluid or susceptibility-current moment kernels |
| [[calculate_heating]] / [[write_energy_record]] | State-isolated corrected diagnostics and versioned companion writer |
| [[zet_in]] | Production-facing signed-k Z interface, selectable backend |
| [[rtsec]] | Existing secant method; not a replacement root algorithm |
| [[fpc(module)]] | Historical all-kinetic FPC distribution/debug routines |

The new numerical modules follow compiler-promoted production real/complex
kinds and strict IEEE behavior where required. Do not assume nrtype::DP is
binary64 under every promoted build. The Faddeeva C ABI explicitly uses the
ISO C double-complex interface. Shared kinetic includes avoid independently
edited copies of the equations; fixed128 compiler-sensitive arithmetic is
protected by exact tests.

<a id="verification"></a>

## Verification

```bash
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

Tests include independent cold fluid/current solutions, NHDS/Stix tensor
references, SciPy/Faddeeva and optional arbitrary precision, componentwise
susceptibility checks, roots, field/nullspace/Faraday checks, energy derivatives,
reader markers, wrapper roundtrips and checked Fortran builds. A passing suite
is evidence for those cases, not universal validation of collisional/FPC or
strong-mode interpretations. Exact expected data should not be edited to hide
numerical changes. [Historical reproduction](reproducibility.md) identifies the
pinned comparisons and settings.

<a id="local-documentation-build"></a>

## Local documentation build

The current FORD configuration points to `docs/`, generated from a fixed staging
list in `.github/workflows/doc.yml`. It is not automatic discovery of every root
Markdown file. The source manual's new pages therefore need additional staging
in a **temporary preview** to inspect the whole hierarchy; this task does not
modify the deployment workflow, page-order configuration or publish anything.

With FORD/Graphviz installed, the existing staged workflow invokes
`ford ford_project.md`. A bare invocation before staging static pages can build
API output while omitting the manual. Current staging also rewrites `.md)` links
in selected files but does not generally resolve anchors, all historical report
links, missing assets, or root-relative example downloads. These build/link
limitations are recorded in the [local manual report](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_5B1_MANUAL_REPORT.md),
with reproducible checks and an expanded temporary preview. Do not claim a new
source page is published merely because its Markdown is committed.

Future infrastructure work should add an explicit staging/navigation policy,
copy required assets once, validate fragment/file/API links, check public input
existence/default tables and run bounded examples. It must remain separately
reviewable from scientific manual content. Deployment remains unchanged here.

<a id="historical-development-evidence"></a>

## Historical development evidence

The following reports describe the state **at their development milestone**.
Statements such as “heating is deferred” or “warm species with a cold reference
are rejected” may be superseded. Use [Models](models.md), [Heating](heating.md)
and [Output](output.md) for the current supported contract. The reports retain
scientific derivations, independent benchmarks and implementation rationale.

| Topic | Historical report / evidence |
|---|---|
| Dual normalization | [Original implementation audit](https://github.com/kgklein/PLUME/blob/main/tests/NORMALIZATION_IMPLEMENTATION_HISTORY.md) |
| Cold kernel convention | [Cold susceptibility design](https://github.com/kgklein/PLUME/blob/main/tests/COLD_SUSCEPTIBILITY_DESIGN.md) |
| Mixed dispersion, warm reference | [Mixed dispersion report](https://github.com/kgklein/PLUME/blob/main/tests/MIXED_DISPERSION_REPORT.md) |
| Fully cold dispersion | [Fully cold report](https://github.com/kgklein/PLUME/blob/main/tests/FULLY_COLD_DISPERSION_REPORT.md) |
| Auxiliary scale mathematics | [Auxiliary audit](https://github.com/kgklein/PLUME/blob/main/tests/AUXILIARY_THERMAL_NORMALIZATION.md) |
| Cold-reference mixed production | [Cold-reference report](https://github.com/kgklein/PLUME/blob/main/tests/COLD_REFERENCE_MIXED_REPORT.md) |
| Adaptive harmonics/cache | [Bessel report](https://github.com/kgklein/PLUME/blob/main/tests/ADAPTIVE_BESSEL_CACHE_REPORT.md) |
| Z accuracy audit | [Plasma dispersion audit](https://github.com/kgklein/PLUME/blob/main/tests/PLASMA_DISPERSION_AUDIT.md) |
| Faddeeva integration | [Faddeeva report](https://github.com/kgklein/PLUME/blob/main/tests/FADDEEVA_BACKEND_REPORT.md) |
| Collisional zz | [Longitudinal correction](https://github.com/kgklein/PLUME/blob/main/tests/COLLISIONAL_ZZ_CORRECTION.md); not general collisional validation |
| Mixed scans | [Scan report](https://github.com/kgklein/PLUME/blob/main/tests/MIXED_SCAN_REPORT.md) |
| Electromagnetic fields | [Field report](https://github.com/kgklein/PLUME/blob/main/tests/MIXED_EIGENFUNCTION_REPORT.md) |
| Species moments | [Cold-fluid/kinetic moments](https://github.com/kgklein/PLUME/blob/main/tests/COLD_FLUID_EIGENFUNCTION_REPORT.md) |
| Heating audit and correction | [Audit](https://github.com/kgklein/PLUME/blob/main/tests/MIXED_HEATING_ENERGY_AUDIT.md), [corrected implementation](https://github.com/kgklein/PLUME/blob/main/tests/MIXED_HEATING_ENERGY_REPORT.md) |
| Documentation/output inventory | [Stage 5A](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_DOCUMENTATION_OUTPUT_AUDIT.md) |
| Targeted reader/workflow fixes | [Stage 5B.0](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_5B0_CORRECTIONS.md) |

Companion audit JSON files remain alongside these reports. None of their frozen
historical measurements is rewritten to match current documentation or tests.
