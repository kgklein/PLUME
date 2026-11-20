# Stage 5B.1: authoritative user manual and scientific reference

## 1. Scope and verified repository state

Task-start/current HEAD: `a853ca0ed4882e9738e4fcc29b8efb1cf4f64013`.
Tracked files were clean at task start. Existing untracked research artifacts,
Stage 5A/5B.0 evidence and previous numerical results were preserved.
This is a **documentation-only scientific update** plus an isolated read-only
manual checker. No numerical Fortran, Python reader behavior, input examples,
output schema, golden expectation, build rule or deployment configuration was
changed. The YAML configuration portion of ford_project.md remains byte-identical;
only its landing-page body changes. Nothing was committed, pushed or deployed.

The published [homepage](https://kgklein.github.io/PLUME/) was checked during
this task and still contained the inherited parallelized/relativistic claim.
This source update corrects that statement; it does not imply that the live site
has changed. Published-site inventory findings from Stage 5A remain historical
evidence, not a new certification of every external page.

The source baseline was established from the current namelists/validation,
composition and normalization dispatch, species tensors, field/moment kernels,
corrected heating snapshots/derivatives, writers and Python readers. The
normalization, cold/mixed, auxiliary, Bessel, Faddeeva, field, moment, heating and
reader reports supplied mathematical/validation context. Superseded statements
about deferred moments/heating/cold references were not carried into the current
capability tables.

## 2. Documentation structure and page inventory

Sixteen manual/appendix pages form the scientific user guide, with a separate
FORD landing-page body. Electromagnetic fields and species moments share one
page; API guidance and historical-report indexing share another. This avoids a
separate tiny page for every development stage or routine.

| Source | Action | Current purpose |
|---|---|---|
| README.md | Rewritten | Current overview, workflow-oriented manual index, verified scope, citation/license |
| INSTALL.md | Rewritten | Actual compiler/runtime/optional Python/test/doc requirements; execution and verification |
| tutorial.md | Rewritten | Bounded three-point mixed-plasma Quick Start, annotated input, actual output and Python |
| input.md | Rewritten | All 89 distinct inputs, types/defaults/units/dependencies/validation/examples, experimental layouts |
| normalization.md | Rewritten | Authoritative physical definitions, anisotropic rho/d conversion, cold references, output units |
| models.md | New | Bi-Maxwellian/cold/mixed models, classification, equilibrium warnings and capability matrix |
| modes.md | New | Direct/supplied/map/serial/2D/repeated/radial workflows and actual root/flag behavior |
| eigenfunctions.md | New | Existing cofactor Ex=1 fields, Faraday factors, kinetic/cold moments, current and failure isolation |
| heating.md | New | Corrected/legacy absorption, all seven channels, complete energy, derivative and applicability definitions |
| output.md | Rewritten | Exact variable-species-count schemas, companion joins/statuses, current reader/wrapper contract |
| examples.md | New | Twelve-lesson index using actual inputs; missing curated lessons explicitly labeled |
| numerics.md | New | Tensors, root/continuation, harmonic bounds/cache, Z and cold-fluid response |
| limitations.md | New | Single restriction/troubleshooting reference; unsupported versus unavailable versus failed |
| reproducibility.md | New | Three independent legacy settings versus an actual historical source revision |
| development.md | New | API/dependency map, verification, build caveats and historical evidence index |
| README-JETPLUME.md | Rewritten | Current all-kinetic FPC appendix; corrected routine/file descriptions and explicit scientific caveats |
| ford_project.md body | Rewritten | Accurate nonrelativistic serial overview; manual as user entry point rather than generic procedures |

Cross-page links point to the central definition rather than repeating full
normalization/heating derivations throughout every table. Explicit source
anchors make fragments work without changing FORD configuration. Relative
manual-home links use the directory index, avoiding the workflow's erroneous
README.html target while remaining usable in the repository.

The former root normalization implementation report is archived verbatim after a
historical-location notice in tests/NORMALIZATION_IMPLEMENTATION_HISTORY.md.
Its original source body matches the task-start file exactly. The developer
history table links existing reports as milestone evidence; it warns that their
restrictions can be superseded. No prior report or audit JSON was rewritten.

## 3. Scientific-content corrections

The guide now distinguishes ordinary species susceptibility chi_s, dielectric
epsilon=I+sum chi_s, electromagnetic Lambda, raw det(Lambda), historical direct
stdout scaling and complex mode frequency. It fixes the old assertion of
unrestricted propagation and separates dispersion support from diagnostics.
The exact-parallel mixed/cold limits are not advertised as all-kinetic support.
Fully cold perpendicular single-wavevector dispersion is separate from the
prohibited kpar=0 scan limit and the Ex=0 ordinary-mode gauge restriction.

Physical drifts remain parallel and Alfvén-normalized, charge ratios signed,
kinetic temperatures bi-Maxwellian and cold selectors exact zero. Neutrality,
current balance and reference drift are warnings only, without automatic density,
frame or drift adjustment. Nonstandard homogeneous equilibria may require
external maintenance.

The guide explicitly retains the current logarithmic map inconsistency, with
its actual origins/increments and lost imaginary signs. Linear maps are
recommended; no numerical correction or compatibility switch is invented.
Map candidate output precedes later filters, and nroot_max is not a guaranteed
.roots row cap or an enabled least-damped sorting policy.

FPC remains separate and restricted. The appendix corrects swapped option6/7
prose and obsolete df1-versus-dfs names, and documents the literal internal
correlation weight without claiming a newly validated dimensional-power
normalization. Its header stores w*sqrt(beta)/Kpar, hence physical phase speed
in v_A units divided by sqrt(alpha_R), not an unconditional Alfvén phase speed.
Fixed FPC distribution harmonics and debug quadrature limitations are explicit.
Some writers already zero/saturate values; reader NaN preservation cannot undo
that loss. No FPC calculation or normalization was changed.

## 4. Complete input-reference coverage

The read-only checker independently extracts the current Fortran namelists:
**89 distinct names in 14 definitions**. input.md has exactly one documented
parameter anchor/record per distinct name, no missing or extra names, organized
as params, species, guess, maps, scan_input, fpc, radial_input and the seven
alternative k_range layouts. The eight control tables plus radial-layout table
all render: **nine input tables**, not raw collapsed Markdown pipes.

Every record includes type, explicit default or “no safe default,” definition,
normalization, allowed/intended values, applicable mode/dependencies, actual
validation and an existing example or explicit lack of a curated one. Optional
status is conditional on the workflow: values with no safe initialization must
be supplied when their group is used. Compiler zeroed storage is not presented
as a default. Incomplete finite/count/storage/Kn/warm-density checks are not
misrepresented as comprehensive enforcement.

The guide documents the restrictive lowercase standalone indexed headers,
required output directories, writable input/scratch directory, short-path
limitation and independent writeOut/use_map semantics. It covers all new
selectors without renaming them or inventing public controls for solver_* state,
Z tolerance, auxiliary beta, nbesmax or diagnostic internals. The wrapper's
supported scalar scope, optional field preservation, dataName handling and
explicit rejection of unsupported controls are described against Stage 5B.0.

## 5. Normalization coverage

The physical species1 reference, positive reference frequency convention,
perpendicular thermal gyroradius, inertial length, reference Alfvén speed,
thermal-speed definitions, beta and signed mass/charge/density/drift ratios
are defined together. Conversion retains sqrt(beta_R*alpha_R), not just
sqrt(beta_R). d_ref temperature is distinguished from rho_ref inverse ratio.
Beta-only scans in different conventions are not claimed physically equivalent.

The unchanged normalization_rho_ref.in/normalization_d_ref.in pair was rerun.
For beta_R=.4, alpha_R=1.5, a=1e-4, the length ratio is sqrt(.6). All three
printed frequencies match after coordinate conversion; the first is
0.06024424 - 0.0004186511 i. This is a verified example result, not a theoretical
frequency invented for the manual.

Cold-reference beta=0 and tau1=0 never introduce a physical thermal reference.
The internal first-warm auxiliary scale, isotropic normalization anisotropy and
cancellation from physical response/output are explained without treating the
scale as a temperature, species or public input. Moment velocities use c*deltaEx/B0,
density uses each population's n0, and E/B use deltaEx in Gaussian units.

## 6. Diagnostic definitions

Fields use the existing complete-wave-matrix cofactor method and Ex=1 phase/
amplitude convention. Faraday's Gaussian c and native-normalization factors are
explicit; an Ex=0 mode can have a valid root but unavailable fields. Kinetic
moments use susceptibility current plus continuity; cold moments solve the
normalized momentum/continuity equations. Current includes density perturbation
advecting equilibrium parallel drift. Species ordering and isolated failures
are preserved.

Corrected heating defines raw Im(E-dagger chi_s E), dimensionless Ps=Qabs/H,
the approximate weak-mode relation to -gamma/omega_r and its limits. All seven
historical channels and ordering are retained, including why four n0 terms close
to full n0, why selected 0/+1/-1 do not generally close to species total, and why
summing all seven double-counts n0. Cold kinetic slots remain NaN.

Complete formal H includes vacuum electric, real-frequency magnetic and all
physical species Hermitian dielectric derivatives. The stationary cold
oscillatory-bulk-energy identity is documented with its zero-drift scope.
Signed fractions and noncancelling f_E_abs have distinct meanings; negative
formal energy is not automatically numerical failure. Real-frequency snapshots,
centered/Richardson matrix-plus-projection convergence, analytic cold derivatives,
state isolation and separate numerical/applicability statuses are described.
Collisional corrected heating and mixed FPC are not newly enabled.

## 7. Output-schema and reader coverage

The output page traces every current family: direct text, five-column maps,
6+6S candidate roots, all scan flag/low_n combinations, component-major tensor
files, the 27+15S numeric-plus-quoted-string energy companion, normalization and
Bessel metadata, gyro/Cartesian/3D FPC, debug eigen/moment files and experimental
radius-last records. It supplies zero-based formulas and one-based physical
species conventions without confusing the two.

It explicitly documents parameter offsets, all seven channel indices, original
species order, real/imaginary pairing, no frequency columns in tensor files,
2D flag rules, legacy low_n=false orphan records and the non-energy/unassigned
historical Ew field. No new energy terms are inserted into old heating blocks.

Energy joins use source_file and scan/outer/internal point/mode coordinates,
not raw row equality with hidden substeps. Failed written roots retain records;
interruption is not atomic across buffered modes and the companion. The current
sweep, contour and energy readers are distinguished from debug/radial/FPC
schemas. API-only scan_index selection, explicit schema flags, cold NaNs and
optional plotting-only FPC NaN conversion match the current Python behavior.
The limitations page also distinguishes the valid contour file-reader schema
from reference-thermal resonance overlays whose helper divides by beta_R and
is undefined for cold references. No overlay extension is advertised.

## 8. Existing examples and executable snippets

**Sixteen unchanged input files** were run in isolated temporary directories.
All processes completed within the bounded run timeout. Fifteen had no failed
scan points; mixed_scan_cold_reference.in had one explicit failed written point
with subsequent recovery. This is documented as a failure/recovery lesson rather
than an all-success introductory example.

| Executed input(s) | Finding used in the manual |
|---|---|
| mixed_fields_warm_reference | Three finite rows; first w=.06153842-5.564326e-6i; 46-column eigen-only schema |
| mixed_fields_cold_reference | Successful bounded cold-reference field/moment scan |
| cold_fields_parallel | First w=.05820986 with roundoff imaginary part; corrects the audit's mismatched example value |
| cold_proton_electron, cold_proton_electron_perpendicular | Supplied single-wavevector cold roots; no nonexistent standalone .roots file promised |
| normalization_rho_ref, normalization_d_ref | Equivalent physical conversion; all three printed frequencies match |
| heating_corrected_kinetic | First H=20342.68004068636; f_E=4.9219490059e-5; status/applicability zero |
| heating_corrected_warm_reference | First H=20384.800070073816; f_E=4.9154412448e-5; cold heat NaN |
| heating_corrected_cold_reference | First H=20516.87918170645; f_E=4.8980504297e-5; beta_ref=0 |
| heating_corrected_cold | First H=20607.629616160644; f_E=4.8875452236e-5; kinetic heat unavailable |
| mixed_moments_multispecies | Finite moments in original four-population order |
| heating_legacy | Explicit historical numerical/diagnostic options; no corrected companion promised |
| mixed_scan_warm_reference, fully_cold_scan | Existing serial examples without failed points in this replay |
| mixed_scan_cold_reference | Explicit NaN point and recovery; exit zero not equated with full success |

The twelve-lesson page labels missing curated anisotropy-instability, whistler/R
and small-2D inputs as **not created**, even though independent developer tests
exist. A separate moderate-a artificial-mass-ratio energy benchmark is not
misattributed to the supplied physical-mass-ratio cold example. Large historical
map/grid and numerical stress inputs remain secondary, not introductory.

**Nineteen distinct downloadable input targets** are linked and verified to
exist in this checkout. External journal/repository HTTP links are not all
claimed live-tested. Input links use repository source URLs because the existing
stager does not bundle input downloads.

All **seven executable documentation checks** pass: six Python snippets from
tutorial/heating/output, plus the complete annotated Quick Start namelist.
The latter produces a **byte-identical mode file** to the supplied bounded
example. Literal placeholder FPC filenames are explicitly comments, not
advertised runnable public examples. No input files were rewritten.

## 9. Local builds, math and presentation inspection

FORD was not initially installed in the normal Python environment. FORD 7.0.11
and its dependencies were installed only under /tmp/plume-5b1-doc-python.
Graphviz and headless Google Chrome were available. Production requirements,
Python environment configuration and deployment files were not changed.

Two builds were inspected in temporary copies:

| Build | Result and limits |
|---|---|
| Existing workflow staging, unchanged build script | FORD succeeds, but ten manual/appendix pages are omitted and 58 manual link/target occurrences remain unresolved (55 internal and three repository URLs) |
| Expanded preview-only staging of root manual pages/reports/assets | FORD succeeds; all 16 manual pages plus landing page present; zero missing manual links/anchors/assets/API targets |

The expanded preview has **35 rendered manual tables**, including all nine
input tables, **82 mathematical script blocks**, and **37 distinct referenced
API/list destinations** that resolve. All 89 parameter anchors survive rendering.
The checker verifies actual tables, not just source pipe notation.

Headless Chrome independently typeset **all 82 expressions across nine scientific
pages**, with no MathJax/merror elements. Screenshots of normalization, heating,
manual home and input tables were inspected. This checks actual MathJax output,
not merely that backslashes survived Markdown. An input-table blank-line problem
found during visual review was corrected and added to the checker's validation.
ASCII multiplication asterisks in reference-table cells are escaped so they
cannot silently turn physical units/formulas into Markdown italics.
The tested browser loads FORD's existing CDN MathJax; offline availability of that
external resource is not claimed.

Explicit manual anchors solve content fragment IDs without configuration changes.
The current template still duplicates metadata/body titles in places, favors
API navigation and lists historical reports prominently. The manual's own home
navigation provides the workflow hierarchy; template/sidebar/staging improvements
remain Stage 5B.3 work, not silently applied here.

Artifacts remain locally under:

```
/tmp/plume-5b1-wkmgita9/workflow-build/docs_out/
/tmp/plume-5b1-wkmgita9/preview-build/docs_out/
/tmp/plume-5b1-wkmgita9/preview-checks.json
/tmp/plume-5b1-wkmgita9/workflow-checks.json
/tmp/plume-5b1-wkmgita9/rendered/results-final.json
/tmp/plume-5b1-wkmgita9/snippet-checks.json
/tmp/plume-5b1-wkmgita9/examples.json
```

Open the expanded `preview-build/docs_out/page/index.html` for the manual,
not an assumption that its contents have been deployed.

## 10. Known build/link failures, deliberately not repaired

The unchanged staging list omits the nine new guide pages and the existing
JET-PLUME appendix. Its simple `.md)` rewriting does not handle `.md#fragment`
links generally and also rewrites external GitHub Markdown targets incorrectly;
it copies only selected images and reports. Its successful
FORD exit is therefore not a complete-manual readiness check. The new read-only
checker correctly exits nonzero on this staged site.

The expanded preview fixes staging **only inside the temporary copy** for
inspection. Across the entire 256-HTML-file preview (including historical/API
pages), a broader scan still finds **27 unresolved occurrences outside the
current manual**:

- 13 FORD-generated interface anchors in poly/assert pages.
- Eight historical report JSON resources not staged beside their HTML.
- Five old relative-path assumptions in the relocated normalization history.
- One old heating-report output-page/fragment link.

Those historical paths are not a claim that current manual links failed; the
manual points to the repository evidence index. The relocated report explicitly
states its original location and links the pinned original. Source-history links
and flattened website layouts need a separate staging policy.

Build warnings also remain: unknown `srd_dir`, `base_url` and project-level
`ordered_subpage` metadata; source-comment metadata parsed from rho_ref/d_ref
text; and an existing alternate-comment parse warning in unlinked
gk_disp_gf.f90. The production build does not link that old file. No numerical
source comment or FORD/deployment setting was changed to silence warnings.
Optional ADS citing-paper generation was skipped without a token in local builds.
No deployment step, GitHub Pages operation or write to the published site ran.

## 11. Reproduction commands

Production and numerical regression checks:

```bash
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

On this machine optional mpmath is under /tmp/plume-z-python, so both full runs
used `PYTHONPATH=/tmp/plume-z-python` as well. Initial baseline: **230 passed in
135.596 s**. Final run after the manual work: **230 passed in 140.858 s**, with
no failures/errors/skips. Checked Fortran builds and exact historical source
comparisons are included; no golden values or existing test expectations changed.

Local source/manual checks can be rerun without deploying:

```bash
python3 tests/check_user_manual.py
python3 tests/check_user_manual.py --site /tmp/<preview>/docs_out --output /tmp/manual-checks.json
```

The checker uses the existing read-only namelist extractor and standard Python
HTML parsing. It certifies local source targets, rendered tables/math markup,
manual links/fragments/assets and linked generated API destinations. It does not
certify every external reference or mathematically validate arbitrary new text.
Browser verification is an additional check.

For a fresh expanded local preview with FORD installed, the following stages only
a temporary copy and does not edit project configuration:

```bash
python3 - <<'PY'
from pathlib import Path
import re, shutil, subprocess, tempfile
from urllib.parse import urlsplit

def local_links(text):
    def convert(match):
        target=match.group(1); url=urlsplit(target)
        if url.scheme or url.netloc: return match.group(0)
        return ']('+re.sub(r'\.md(?=[#)]|$)','.html',target)+')'
    return re.sub(r'\]\(([^)]+)\)',convert,text)

root=Path.cwd()
preview=Path(tempfile.mkdtemp(prefix='plume-manual-preview-'))
shutil.copytree(root/'src',preview/'src')
shutil.copy2(root/'ford_project.md',preview/'ford_project.md')
docs=preview/'docs'; docs.mkdir()
for source in [*root.glob('*.md'),*(root/'tests').glob('*.md')]:
    if source.name=='ford_project.md': continue
    title=re.search(r'^# (.+)$',source.read_text(),re.M)
    name='index.md' if source.name=='README.md' else source.name
    text=local_links(source.read_text())
    (docs/name).write_text('title: '+(title.group(1) if title else source.stem)+'\n\n'+text)
for name in ['PLUME_logo.png','qrcode_plume_github.png','Jet-Plume_Logo.svg']:
    shutil.copy2(root/name,docs/name)
subprocess.run(['ford','ford_project.md'],cwd=preview,check=True)
subprocess.run(['python3',str(root/'tests/check_user_manual.py'),'--site',str(preview/'docs_out')],check=True)
print(preview/'docs_out/page/index.html')
PY
```

To reproduce the **existing workflow** instead, execute its build-only staging
shell block in a separate copy with FORD on PATH; do not execute its deployment
step or run its source-rewriting sed commands against a research checkout.

## 12. Remaining work and implementation boundaries

1. **Stage 5B.2:** create bounded missing curriculum inputs (instability,
   whistler, small 2D, moderate-a energy), independently identify branches,
   test expected outputs/reader schemas and automate example/default coverage.
   Do not promote a stress scan just because its process exits zero.
2. **Stage 5B.3:** make staging/navigation/assets/report JSON/source-relative
   link policy explicit; handle local fragments without rewriting external source URLs, API anchor issues, documentation
   warnings and optional citation generation; run successful local preview/link
   checks in CI before any separately authorized publication.
3. **Stage 5B.4:** cross-audience scientific review, tutorial usability checks,
   final rendered copy-editing and documentation/example maintenance procedures.
   Keep legacy source/settings provenance and availability distinctions intact.
4. All-kinetic parallel physics, perpendicular warm limits, map-equation changes,
   collisional off-diagonal correction, dimensional FPC normalization/closure and
   mixed FPC remain separate scientific tasks. Documentation does not implement
   them or infer correctness from unrelated passing tests.

## 13. Complete changed-file inventory and preservation

**Modified existing documentation (8):** README.md, INSTALL.md, tutorial.md,
input.md, normalization.md, output.md, README-JETPLUME.md, and the body of
ford_project.md.

**New manual pages (9):** models.md, modes.md, eigenfunctions.md, heating.md,
examples.md, numerics.md, limitations.md, reproducibility.md, development.md.

**New documentation/audit support (3):**

- tests/NORMALIZATION_IMPLEMENTATION_HISTORY.md — archived old normalization report.
- tests/check_user_manual.py — isolated read-only source/site verification.
- tests/PLUME_5B1_MANUAL_REPORT.md — this report.

Final protected-file comparison covers 228 pre-existing tracked/test source files:
220 unchanged, and exactly the eight intended existing documentation edits.
Production Fortran, Python reader/wrapper code, input examples, build/deployment
configuration, old reports/measurements and test expectations remain unchanged.
The ford_project configuration prefix and archived normalization source were
checked against the original bytes. Existing assets are preserved. HEAD remains
the task-start commit; all work is local and uncommitted.
