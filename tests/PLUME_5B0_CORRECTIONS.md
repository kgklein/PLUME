# Stage 5B.0: targeted reader and workflow corrections

## Summary and scope

Task-start commit: `d012b9c480f80745145b8dbf6d0bfbf89dbd110e`.
Tracked files were clean at task start; unrelated untracked research files and
the Stage 5A inventories were retained. This stage corrects Python interpretation
and serialization, and adds input/preflight warnings. It changes no susceptibility,
normalization, Z backend, harmonic selection/cache, root iteration, heating,
FPC numerical method, map coordinate equation or numeric output format.
Everything remains local and uncommitted.

The confirmed fixes are:

- Contour heating detection now uses the selected input schema and flags.
  Eigen-only species parameters are not heating.
- Gyro/Cartesian FPC readers retain NaNs, infinities, signed finite samples,
  physical zeros and shapes. NaN-to-zero conversion requires an explicit
  visualization option.
- Wrapper `writeOut` and `use_map` are independent. The analogous accidental
  coupling of `positive_roots` to `loggridw` is also removed.
- Supported scalar namelist controls survive wrapper roundtrips, including
  all omitted map, FPC and tensor fields identified by Stage 5A.
- Unsupported all-kinetic parallel and ignored map diagnostic requests warn.
- Logarithmic maps are characterized as an internally inconsistent historical
  implementation, not redefined as a legitimate physical-bound convention.

## Logarithmic map evidence and compatibility decision (before changes)

The current map equations mix two conventions. For real-axis node i=0..nr,

    dr = (log10(omf)-log10(omi))/nr
    wr_i = 10**(omi+i*dr)

Thus wr_0=10**omi and wr_nr=10**omi*(omf/omi), not omi and omf.
For imaginary node j=0..ni,

    di = (log10(abs(gamf))-log10(abs(gami)))/ni
    wi_j = -10**(abs(gami)+j*di)

Both imaginary input signs are discarded; every log imaginary node is negative.
The final/initial magnitude ratio is abs(gamf/gami), but the starting magnitude
is 10**abs(gami), rather than abs(gami). The audit's endpoint examples follow
exactly. Negative real endpoints and zero imaginary magnitudes are not valid
arguments to the logarithms.

Neither physical-frequency bounds nor exponent bounds explain both the origins
and increments consistently. An exponent convention would instead need
(omf-omi)/nr and an explicitly defined imaginary exponent/sign convention.
Current public descriptions call these frequency bounds, and the stock map
examples use linear grids. History establishes how the mismatch arose:
`86ddd30` (2023-01-18) used exponent origins **and linear increments**;
`ae14a90` (2025-01-08) changed the increments to logarithms of bound ratios,
retained the real exponent origin, and changed the imaginary origin/sign to
`-10**abs(gami)`. That numerical change appears in a commit whose message
says no functionality changed. The present inconsistency is historical but
is not a coherent intentional exponent convention. We cannot infer the
author's intended physical convention from the commit message.

**Decision for this maintenance stage:** preserve all map equations, coordinates,
root methods, and numeric file schemas. Add characterization tests and an input
warning explaining the historical convention. A future, separately approved
correction should define physical positive real endpoints and same-sign,
nonzero imaginary endpoints, interpolate logarithms of magnitudes, and preserve
an explicit imaginary sign. It would change sampled determinants, selected
roots, maps, and continuations; exact historical log-grid reproduction therefore
requires an explicit compatibility policy or the old source revision. Existing
legacy Z/Bessel/heating switches do not select a map convention.

### Numerical probes

The production map writer returns five columns: real frequency, imaginary
frequency, log10 determinant magnitude, real determinant and imaginary determinant.
`nr=ni=2` gives 3x3 nodes; the blank block separators are unchanged.

| Requested bounds | Linear nodes | Historical logarithmic nodes |
|---|---|---|
| wr=.1..1 | .1, .55, 1 | 1.2589254, 3.9810717, 12.589254 |
| gamma=-1e-4..-1e-2 | -.0001, -.00505, -.01 | -1.0002303, -10.002303, -100.02303 |
| gamma=+1e-4..+1e-2 | +.0001, +.00505, +.01 | same negative nodes as the previous row |

The real start is too high by `10**omi/omi`; the imaginary starting magnitude
is too high by `10**abs(gami)/abs(gami)`. Grid ratios are correct only relative
to these displaced origins. Tests exercise each log flag independently and
together, both imaginary signs, and the checked executable. The decimal output
is compared with independent formula evaluation at `rtol=6e-7`, matching
the existing six-place map-coordinate precision. Checked/production printed
coordinates agree exactly; determinant columns agree at `rtol=1e-12`.

All bundled public map examples use linear flags. This avoids changing their
established behavior, but does not prove that no private legacy workflow uses
the inconsistent logarithmic equations. The characterization tests intentionally
pin current behavior and must be reviewed explicitly in any future map fix.

## Unsupported-mode compatibility decision (before changes)

Warn rather than reject the previously accepted exact-parallel all-kinetic and
ignored map diagnostic requests. This preserves legacy accepted-input/exit
behavior and numerical results. Warnings occur during input/scan preparation,
not each dispersion/root evaluation. No parallel kernel, wavevector floor,
map diagnostics or new output columns will be introduced. Existing hard
restrictions (for example corrected map heating) remain in force.

## Reader schemas

The contour reader accepts serial `option=1` and two-dimensional `option=2`
mode files with explicit matching input metadata. It does not guess a schema
from surplus numeric columns. Let S be species count, E the eigen flag, P the
heating flag, and L the low_n flag:

| Block | Zero-based start | Length |
|---|---|---|
| Plasma/wavevector/frequency prefix | 0 | 6 |
| Fields plus species moments, if E | 6 | 12+8S |
| Species totals, if P | F | S |
| Seven channels, if P and L | F+S | 7S |
| Species parameters | F+K | 6S |
| Historical Ew slot, only E and P and not L | F+K+6S | 1 |

Here F=18+8S if E, otherwise 6. K=8S if P and L, S if P and not L,
otherwise zero. Parameters are six values per species in original input order.
The exact expected width is F+K+6S plus the conditional old Ew slot. Ew is not
interpreted as corrected wave energy; use the existing `.energy` companion.

The seven channels use species-major stride seven, in order: TTD yy, TTD yz,
LD zy, LD zz, full n=0, transverse n=+1, transverse n=-1. Heating-only starts
at column seven (one-based). Eigen-only heat_start is None and no P_species
keys are returned. NaNs remain NaNs in all blocks.

Serial inputs with unequal diagnostic flags require `scan_index=N` (one-based)
in `read_dispersion_output`. Two-dimensional output uses scan_input_2 flags.
Serial namelists are sorted numerically rather than lexicographically.
Incorrect widths, unknown file families and ambiguous layouts raise useful
ValueErrors. Ragged historical low_n=false orphan records are rejected
explicitly; rectangular low_n=false records use the specified schema. No
old writer or column position changes. The generic `load_plume_sweep` and
`load_plume_energy` implementations are unchanged.

The resonance plotter CLI does not yet expose the API's scan_index selector.
For heterogeneous serial diagnostic flags, use the explicit reader API or a
matching single-scan input; the CLI refuses ambiguity instead of plotting a
plausible but wrong layout.

### FPC missing samples

`loadlinfpccepar`, `loadlinfpcceperp`, `loadlinfpccart`, `loadlinfpc3d` and their
distribution-loading helpers accept the keyword-only option
`visualization_nan_to_zero=False`. Default conversion is ordinary float parsing.
The optional mode replaces only NaNs; it preserves infinities and all finite
samples. The 3D reader's old finite tail/saturation clipping is also removed
from its default loading. The private block reader has preserving defaults.
Velocity axes, projection order, distribution units, resonance metadata,
species naming heuristics and physical FPC calculations are untouched.

Malformed-file retry/offset heuristics and hard-coded species labels in these
historical readers are outside this correction. They still need future review;
this stage does not advertise general mixed-plasma FPC support. Existing file
handle ResourceWarnings can appear under unittest; they are not new numerical
failures and broad reader resource/error refactoring is deferred.

## Wrapper scalar serialization

The old line/whitelist loader silently omitted controls. A small shared scalar
schema now covers every field in the six groups the wrapper supports:

| Group | Controls | Examples newly retained |
|---|---:|---|
| params | 23 | low_n, dataName |
| species_N | 6 | all original mass/charge/density/temperature/drift values |
| maps | 9 | loggridg, nr, ni |
| scan_input_N | 10 | tensor |
| guess_N | 2 | both signed frequency components |
| fpc | 12 | vxmin/vxmax, vymin/vymax, vzmin/vzmax |

Cross-checking the Stage 5A parser inventory finds no missing controls within
these groups. Indexed groups retain numeric order, must be contiguous, and are
replaced on repeated loads. Explicitly supplied values always override defaults.
Only omitted controls receive the few defaults required by the wrapper; omitted
Kn stays None and is omitted on write, retaining Fortran's default.

The parser handles comments, inline scalar assignments, D exponents, quoted
strings with spaces/slashes/exclamation marks and doubled quotes, logical
spellings T/F and dotted variants, and case-insensitive names. Existing
dictionaries retain their canonical Fortran logical string representation.
It is deliberately not a general array/expression namelist parser.
Unsupported groups (including radial_input and k_range), unknown controls,
noncontiguous/duplicate indexed groups and unsupported syntax raise ValueError.
Writers validate control names before replacing an existing input file.

An explicit loaded dataName is now authoritative. To redirect a loaded input,
set obj.dataname afterward. write_input's explicit outputname argument still
selects the output name. This differs from the previous loader's deliberate
discarding of dataName and is necessary for faithful roundtrips. Strings are
quoted exactly once on write; repeat loads do not accumulate species/scans.

`set_maps` keeps its original first six positional arguments and adds trailing
loggridg=False, nr=128, ni=128 options. `make_scan` adds trailing tensor=False.
Neither requires existing callers to add arguments. The independent
positive_roots correction is tested across both values of loggridw.

## Warnings and compatibility

Fortran warnings are confined to functions/scan_support, avoiding edits to
compiler-sensitive disprels arithmetic. Equivalent Python warnings occur when
wrapper input checks/write are requested.

| Request | Current action |
|---|---|
| All-kinetic kperp=0 | Warn once per executable run, including scan preflight; retain historical evaluator/result behavior |
| Cold/mixed kperp=0 | No new warning; existing analytic supported paths retained |
| option4 eigen/tensor | Warn once per requesting scan during input validation; flags remain ignored |
| option4 legacy heating | Warn once per requesting scan; heating remains absent |
| option4 corrected heating | Existing hard rejection retained |
| loggridw or loggridg | One map-input warning describing the displaced origins and lost imaginary sign |

Warnings are not emitted in every frequency/root iteration and remain visible
when writeOut is false. No new rejection of previously accepted supported
Fortran inputs is introduced. Text-only stdout comparisons must account for
the new warnings. Existing numeric files and output schemas remain unchanged.

## Validation

The pre-change checkout passed all **209 tests** in **133.239 s**, with optional
high-precision references available. The final clean build succeeded and the
full expanded suite passed **230 tests in 145.073 s**, with **no failures,
errors or skips**. All 21 new tests passed, including the checked builds.
Both updated defect-characterizations and all unchanged numerical tests passed.

| Check | Result |
|---|---|
| Untouched task-start build and baseline suite | 209 passed |
| Clean production build, make clean then make -j2 | Passed |
| Final full unittest discovery | 230 passed, 145.073 s |
| New workflow/reader tests within full suite | 21 passed |
| Checked executable + input warning probe | Passed with bounds checking and IEEE traps |
| Pinned pre-normalization numerical records | Exact equality passed |
| Pinned pre-Stage-4A kinetic scan/eigenfunction files | Exact equality passed |
| Task-start output-file comparisons | 77/77 byte-identical |
| Python compilation and git diff --check | Passed |

Exact tests include
`test_normalization.NormalizationTests.test_existing_inputs_match_unmodified_solver`
against `c6c7616e1fc58eb6118e1f9423a3f0df6d52ef20`, and
`test_mixed_scans.MixedScanTests.test_exact_kinetic_scan_compatibility` plus
`test_mixed_eigenfunctions.MixedEigenfunctionTests.test_exact_legacy_eigen_scan_files`
against `8331925c7953dcc68f29f14c1d7b78dd9198bd5a`.

New `test_workflow_corrections.py` covers:

- 24 independently constructed schemas (S=1,2,4; all E/P/L combinations),
  unique markers, two rows each, NaNs, all seven ordered channels and parameter
  offsets; ambiguous/mismatched/ragged files and second-axis 2D flags.
- All four writeOut/use_map combinations and all independent map flag pairs.
- Roundtrips of kinetic, warm-reference mixed, cold-reference mixed,
  multi-species, map, scan, gyro/Cartesian FPC inputs and six existing example
  fixtures. Parsed scalar values, not whitespace/formatting, are compared.
- Quoted paths, D exponents, logical variants, repeated loads and rejection of
  unsupported controls before an existing file can be overwritten.
- Synthetic gyro, Cartesian projection and Cartesian 3D arrays containing
  positive/negative values, physical zeros, NaNs/infinities and 1e-120/1e120.
- Twelve real diagnostic scan configurations (four compositions times
  eigen-only/heating-only/both), including unavailable cold heat and unchanged
  `.energy` reading; real gyro and Cartesian FPC output with injected NaN/zero.
- Actual log maps, ignored map requests, parallel warning classification, and
  a full checked executable plus parse-only checked warning probe.

The two Stage 5A defect-characterization assertions are converted to the
requested correctness assertions (independent writeOut and preserved FPC NaN).
No numerical golden data or historical exact-regression assertion is changed.
The Stage 5A Markdown/JSON/site inventories remain historical evidence.

### Numerical preservation evidence

In addition to the pinned historical tests, 13 representative runs were compared
with the untouched task-start executable built in an isolated copy. All **77
generated files are byte-identical**:

| Case | Identical files |
|---|---:|
| Kinetic fixed128 / adaptive scans | 4 / 4 |
| Warm-reference mixed fixed128 / adaptive scans | 5 / 5 |
| Cold-reference mixed fixed128 / adaptive scans | 5 / 5 |
| Fully cold fixed128 / adaptive scans | 5 / 5 |
| Historical logarithmic frequency map | 4 |
| Repeated linear maps with ignored diagnostic flags | 6 |
| Legacy finite-collision-rate scan | 3 |
| Gyro / Cartesian FPC calculation | 9 / 17 |

The files include original numeric columns and existing companions/metadata.
Output equality is exact, not a loosened tolerance comparison. The same
comparison remains valid after the clean rebuild because numerical sources
and build flags were not edited.

### Reproduction commands and toolchain

```
make clean
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
PYTHONPATH=tests MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest test_workflow_corrections -v
```

Environment: GNU Fortran/C 11.4.0, GNU make 4.3, Python 3.10.12, NumPy 2.2.6,
SciPy 1.15.3, Matplotlib 3.10.8. Normal Makefile GNU optimization and promoted
double precision are unchanged. Checked tests compile in a TemporaryDirectory
with `-O0 -g -fdefault-real-8 -fcheck=all -ffpe-trap=invalid,zero,overflow`.
The parse-only probe checks unsupported all-warm parallel input without
deliberately executing its known singular kernel under floating-point traps.
Checked log-map calculations use the supported fully cold evaluator.

On this machine optional mpmath resides in `/tmp/plume-z-python`, so full-suite
commands were prefixed with `PYTHONPATH=/tmp/plume-z-python`. No package or
runtime dependency was added. All generated solver output is isolated under
temporary directories. Task logs/baseline comparison artifacts are under
`/tmp/plume-5b0-oilacf80` and `/tmp/plume-5b0-*.log`.

## Remaining issues and documentation-refresh sequence

1. Decide an explicit map compatibility policy before correcting log equations;
   update the characterization tests only with that separately authorized change.
2. Refresh the user input/output references from Stage 5A's inventories, retaining
   the explicit E/P/L schema and distinction between unavailable and zero data.
3. Document generic low_n=false/ragged legacy parsing limits, the API-only
   scan_index selection, old Ew semantics and debug reader restrictions.
4. Curate executable examples and generated schema/default checks before the
   website restructuring/deployment refresh. This stage changes neither example
   files nor site organization/deployment.
5. Leave all-kinetic exact-parallel physics, map diagnostic formats, mixed FPC,
   radial wrapper support, old species-label heuristics and broad input validation
   to separately scoped work. A warning does not make a calculation supported.

## Complete changed-file inventory

| File | Change |
|---|---|
| contours/plot_plume_resonances.py | Explicit scan diagnostic schema, flags/width checks and numeric scan ordering |
| linfpclib/_namelist.py (new) | Supported scalar group/control schema, parser and serializer |
| linfpclib/linfpc.py | Independent switches, lossless supported scalar roundtrips, missing-value readers and targeted warnings |
| src/functions.f90 | Input warnings only: logarithmic maps and unsupported all-warm parallel |
| src/scan_support.f90 | Once-per-run parallel warning and once-per-scan ignored map diagnostics warnings |
| tests/test_documentation_output_audit.py | Replace two obsolete defect-characterizations with corrected behavior assertions |
| tests/test_workflow_corrections.py (new) | Independent markers, roundtrips, actual files, map coordinates and checked-build tests |
| tests/workflow_input_probe.f90 (new) | Isolated parse-only checked warning harness |
| input.md | Exact retained map equations, independent controls, wrapper roundtrip limits and warnings |
| output.md | Explicit reader schema/availability and FPC missing-value behavior |
| README-JETPLUME.md | Targeted FPC missing-value visualization note |
| tests/PLUME_5B0_CORRECTIONS.md (new) | This report |

Final snapshot check: of 224 protected pre-existing files, 215 are unchanged;
eight have the intended edits listed above. A separate `ford_project.md` edit
appeared during the task (rewording its landing-page capability sentence).
It was not made, reverted or included in this task's changed-file inventory.
That unrelated documentation work is preserved. The four new files are listed
explicitly in the table.

The collisional correction, Stage 5A inventories, production numerical kernels,
Makefile, test goldens and input examples remain unchanged. HEAD is still the
task-start commit. No commit or push is performed.
