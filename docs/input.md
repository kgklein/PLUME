<a id="input-reference"></a>

# Input reference

[Manual home](index.md) · [Normalization](normalization.md) · [Modes](modes.md)

This reference covers the **89 distinct case-insensitive input names** used by PLUME.
They are invoked across 14 namelist definitions, 
where seven alternative `k_range` layouts share names.
The tables describe the parsers and checks.
The source of truth is `functions.f90` namelists and`vars.f90` initialization. 

<a id="syntax-and-required-values"></a>

## Syntax and required values

An input file uses `&params ... /`, 
then indexed `&species_1` through `&species_nspec` namelists,
plus maps/guesses/scans/FPC/radial groups as selected. 
Put each indexed header **alone on a lowercase line**, 
without a trailing comment. 
The custom indexed extractor is case-sensitive even though scalar field names are not. Ordinary commas or newlines can separate scalar assignments. 
Quote character values and use `.true.`/`.false.` for clarity. 
Supply all required species values, including a placeholder `alphS` for cold species. 
See the complete working annotated
[Quick Start](tutorial.md#2-read-the-input-in-physical-terms).

In the below table, “No safe default” means there is no explicit initialization establishing a public value: **supply it whenever the group/workflow uses it**. 
A default removes the requirement to write that particular value.

Common public links: [physical parameters](normalization.md),
[species classification](models.md), [scan semantics](modes.md),
[heating convention](heating.md), [numerical backends](numerics.md),
[FPC restrictions](README-JETPLUME.md), [output schemas](output.md).

<a id="global-controls-params"></a>

## Global controls: &params

This group is required for every run. 
Species type and available diagnostics depend on the selected workflow and geometry.
All quantities must be interpreted in the selected normalization.
The first species is the reference.

| Name / type | Default | Definition, units and allowed values | Existing example |
|---|---|---|---|
| <a id="param-heating_convention"></a>`heating_convention` (character(16)) | `corrected` | Diagnostic convention; corrected complete real-frequency energy and initialized same-frequency harmonics; legacy original calc_eigen arithmetic. Units: dimensionless label. Values: 'corrected'  or  'legacy'. | [heating_corrected_cold.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_cold.in) |
| <a id="param-z_backend"></a>`z_backend` (character(16)) | `faddeeva` | Plasma dispersion Z backend, Z=i sqrt(pi) wofz; signed kpar convention retained. Units: dimensionless label. Values: 'faddeeva'  or  'legacy'. | [heating_corrected_kinetic.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_kinetic.in) |
| <a id="param-bessel_mode"></a>`bessel_mode` (character(16)) | `adaptive` | Species-specific warm harmonic policy. Units: maximum absolute signed order N, sum -N..N. Values: 'adaptive'  or  'fixed'  or  'fixed128'. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-bessel_tol"></a>`bessel_tol` (real) | `1e-10` | Ordinary and n^2-weighted Bessel tail tolerance, not guaranteed tensor/root error. Units: dimensionless. Values: finite 0 &lt; value &lt; 1.| [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-bessel_nmax"></a>`bessel_nmax` (integer) | `512` | Adaptive maximum harmonic cap; warning and continue on unresolved tail. Units: absolute order. Values: 1..46338.| [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-bessel_nfixed"></a>`bessel_nfixed` (integer) | `128` | Fixed maximum absolute harmonic order. Units: absolute order. Values: 1..46338.  | See group workflow / tests |
| <a id="param-bessel_diagnostics"></a>`bessel_diagnostics` (logical) | `.false.` | Enable .bessel cache/cutoff/convergence metadata. Units: logical. Values: .true.  or  .false.. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-betap"></a>`betap` (real) | **No safe default** | 8\*pi\*nR\*Tparallel,R/B0^2. Units: parallel reference beta. Values: positive for warm reference; exactly zero for cold reference. `d_ref` requires tauS(1)=betap within 1e-8 relative|[adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-kperp"></a>`kperp` (real) | **No safe default** | Perpendicular wavevector component, k=(kperp,0,kpar). Units: k\*rho_R or k\*d_R. Values: finite; mixed/cold can be exactly zero. Cold/mixed exact parallel available| [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-kpar"></a>`kpar` (real) | **No safe default** | Parallel wavevector component. Units: k\*rho_R or k\*d_R. Values: finite; nonzero in every wavevector scan and warm susceptibility. |  [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-vtp"></a>`vtp` (real) | **No safe default** | rho_ref: reference parallel thermal speed/c; d_ref: reference Alfven speed/c. Units: v_tparallel,R/c or vA,R/c. Values: positive. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-normalization"></a>`normalization` (character(16)) | `rho_ref` | Physical input length/temperature/velocity convention. Units: rho_ref uses perpendicular reference thermal gyroradius; d_ref uses inertial length. Values: 'rho_ref'  or  'd_ref'. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-nspec"></a>`nspec` (integer) | **No safe default** | Number of populations, in original indexed order. Units: count. Values: must be positive; parser lacks explicit complete bounds check. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-nscan"></a>`nscan` (integer) | `0` | Number of serial scans or 2D coordinate axes. Units: count. Values: nonnegative; 2 for option=2; option=4 one scan. |  [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-option"></a>`option` (integer) | **No safe default** | Executable workflow selector. Units: enumeration. Values:  -1 direct; 0 root map; 1 supplied/map roots + serial scans; 2 grid; 4 repeated maps; 5 radial; 6 gyro FPC; 7 Cartesian FPC. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-nroot_max"></a>`nroot_max` (integer) | **No safe default** | Number of supplied guesses / continued branches, not a guaranteed cap on .roots file lines. Units: count. Values: positive; effective storage ceiling numroots=500, not fully input-validated. |[adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-use_map"></a>`use_map` (logical) | **No safe default** | Choose frequency map search versus supplied guesses. Units: logical. Values: .true.  or  .false.. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-low_n"></a>`low_n` (logical) | `.true.` | Select extra n=0,+1,-1 harmonic outputs in historical mode blocks. Units: logical. Values: .true.  or  .false.. | See group workflow / tests |
| <a id="param-collision_type"></a>`collision_type` (integer) | `0` | Collision operator selection. Units: enumeration. Values: 0 collisionless; 1 nonconservative Krook; other values STOP. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-kn"></a>`Kn` (real) | `1e+30` | Neutral mean free path / selected reference length. Units: Kn_d=Kn_rho\*sqrt(beta_R\*alpha_R). Values: positive physical mean free path ratio; incomplete finite/positivity validation. | See group workflow / tests |
| <a id="param-writeout"></a>`writeOut` (logical) | **No safe default** | Verbosity control for screen output, not a prohibition on files. Units: logical. Values: .true.  or  .false.. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-dataname"></a>`dataName` (character(100)) | **No safe default** | Subdirectory under data/ for output. Units: path. Values: existing relative data subdirectory. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-outputname"></a>`outputName` (character(100)) | **No safe default** | Common prefix for mode/maps/metadata outputs. Units: filename component. Values: output basename. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |

<a id="species-species_n"></a>

## Species: &species_N

Supply six fields for every indexed population. `tauS=0` selects cold only in d_ref; any positive value remains kinetic. Species ordering is never changed automatically.

| Name / type | Definition with `rho_ref` | Definition with `d_ref` | Application notes | Existing example |
|---|---|---|---|---|
| <a id="param-taus"></a>`tauS` (real) | Tparallel,R/Tparallel,s. Units: inverse temperature ratio. Values: &gt;0. | 2\*Tparallel,s/(mR\*vA^2). Units: Alfvenic temperature. Values: &gt;=0; exact zero selects cold. | Every species. Exact zero is model selector only in d_ref; negative invalid; no small-temperature threshold. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-mus"></a>`muS` (real) | Reference/species mass ratio mR/ms. Units: dimensionless. Values: positive. | Reference/species mass ratio mR/ms. Units: dimensionless. Values: positive. | Every species. Physical reference value should be1; parser warns only for mu1&gt;1+1e-7, not symmetric/terminal.  | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-alphs"></a>`alphS` (real) | Species perpendicular/parallel temperature ratio. Units: dimensionless. Values: &gt;0. | Species perpendicular/parallel temperature ratio. Units: dimensionless. Values: warm &gt;0; cold ignored. | Every species; cold placeholder irrelevant. Physical alpha_s distinct from auxiliary alpha\*=1; rho_R uses alpha_R.  | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-qs"></a>`Qs` (real) | Reference/species charge ratio qR/qs. Units: dimensionless. Values: nonzero, signed. | Reference/species charge ratio qR/qs. Units: dimensionless. Values: nonzero, signed. | Every species. Gyrofrequency mu/Q signed; plasma prefactor uses Q^2. No neutrality hard constraint.  | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-ds"></a>`Ds` (real) | Species/reference density ratio ns/nR. Units: dimensionless. Values: nonnegative physical density; warm-only validation incomplete. | Species/reference density ratio ns/nR. Units: dimensionless. Values: nonnegative physical density; cold references/cold kernels enforce; warm-only validation incomplete. | Every species. Zero density has unavailable relative-density moments; not physically defined zero perturbation. | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |
| <a id="param-vvs"></a>`vvS` (real) | Equilibrium Us/vA,R along B0. Units: Alfven velocity. Values: finite signed parallel drift. | Equilibrium Us/vA,R along B0. Units: Alfven velocity. Values: finite signed parallel drift. | Every species, including species1 with warning. No forced zero-current frame or automatic charge/current balancing.  | [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in) |

<a id="supplied-frequencies-guess_n"></a>

## Supplied frequencies: &guess_N

Provide one group per requested guess when `use_map=.false.` in root workflows, and for direct option=-1. Only the first frequency is directly evaluated by test_disp.

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-g_om"></a>`g_om` (real) | **No safe default** | Real part of supplied frequency. Units: Omega_R. Values: finite supplied frequency guess. | &guess_j, use_map=false; direct option-1 uses first. Secant roots are refined except direct evaluation. | See group workflow / tests |
| <a id="param-g_gam"></a>`g_gam` (real) | **No safe default** | Imaginary part of supplied frequency; growth positive. Units: Omega_R. Values: finite, signed. | &guess_j. Strongly damped/growing heating can be qualified/unavailable. | See group workflow / tests |

<a id="frequency-maps-maps"></a>

## Frequency maps: &maps

Required for option=0/4 or any map-initialized workflow. Use both log flags explicitly; linear maps are recommended. These flags concern the frequency plane, not wavevector scan interpolation.

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-loggridw"></a>`loggridw` (logical) | **No safe default** | Select the historical inconsistent logarithmic real-frequency map formula; linear recommended. Units: logical. Values: .true.  or  .false.. | &maps, option0/maps. See Modes for the exact origin/increment equations; bounds are not physical endpoints when true. | See group workflow / tests |
| <a id="param-loggridg"></a>`loggridg` (logical) | **No safe default** | Select the historical inconsistent logarithmic imaginary-frequency map formula; linear recommended. Units: logical. Values: .true.  or  .false.. | &maps. Imaginary signs discarded and all resulting nodes negative; see Modes. | See group workflow / tests |
| <a id="param-omi"></a>`omi` (real) | **No safe default** | Initial real-frequency coordinate for a linear map; historical log origin differs. Units: Omega_R. Values: finite; &gt;0 if log. | &maps. Linear endpoints work; log origin currently exponentiates omi rather than log10(omi). | See group workflow / tests |
| <a id="param-omf"></a>`omf` (real) | **No safe default** | Final real-frequency coordinate for a linear map; historical log endpoint differs. Units: Omega_R. Values: finite; &gt;0 if log. | &maps. Bounds ordering/grid positivity not comprehensively validated. | See group workflow / tests |
| <a id="param-gami"></a>`gami` (real) | **No safe default** | Initial signed imaginary-frequency coordinate for a linear map; historical log origin differs. Units: Omega_R, positive growth. Values: finite; nonzero magnitude if log. | &maps. Log mode always negative; intended bounds currently not reproduced. | See group workflow / tests |
| <a id="param-gamf"></a>`gamf` (real) | **No safe default** | Final signed imaginary-frequency coordinate for a linear map; historical log endpoint differs. Units: Omega_R. Values: finite; nonzero magnitude if log. | &maps. Same log-origin issue as gami. | See group workflow / tests |
| <a id="param-positive_roots"></a>`positive_roots` (logical) | `.false.` | Select candidates with real frequency &gt; -1e-10 after map file writing. Units: logical. Values: .true.  or  .false.. | Map-selected continuation roots. Does not mean positive growth; .roots written before final filter; hard gamma&lt;1 filter also exists.| See group workflow / tests |
| <a id="param-nr"></a>`nr` (integer) | `128` | Real-frequency map intervals, generating nr+1 nodes. Units: interval count. Values: positive intended, not fully validated. | &maps. Default128 -&gt;129 nodes; grid has (nr+1)\*(ni+1) rows.  | See group workflow / tests |
| <a id="param-ni"></a>`ni` (integer) | `128` | Imaginary-frequency map intervals, generating ni+1 nodes. Units: interval count. Values: positive intended, not fully validated. | &maps. Default128 -&gt;129 nodes. | See group workflow / tests |

<a id="scans-scan_input_n"></a>

## Scans: &scan_input_N

Supply nscan indexed groups for serial scans or exactly two axes for option=2. Defaults for diagnostics reset false for each group; range/style/count controls have no safe default. See [scan mappings](modes.md#serial-wavevector-scans) for the overloaded swi/swf meanings.

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-scan_type"></a>`scan_type` (integer) | **No safe default** | Select coordinate inside scan style. Units: enumeration. Values: style-1:0 vector,1 angle,2 magnitude; style0:0 kp,1 kz,2 beta,3 vtp; species:0 tau,1 mu,2 alpha,3 Q,4 D,5 U. | &scan_input_j. Vector swi/swf are final kp/kz; angle uses swi as terminal angle; magnitude uses swf; not generic two endpoints. | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-swi"></a>`swi` (real) | **No safe default** | Scan range_i; compound vector/angle overload. Units: native coordinate; angle degrees. Values: coordinate-dependent valid initial/target value. | &scan_input_j. For style-1 vector: target kperp; angle: target theta; initial physical point comes from current plasma state.  | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-swf"></a>`swf` (real) | **No safe default** | Scan range_f; compound vector/magnitude overload. Units: native coordinate; angle degrees. Values: coordinate-dependent valid final/target value. | &scan_input_j. For style-1 vector: target kpar; magnitude: terminal signed  or k or .  | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-swlog"></a>`swlog` (logical) | **No safe default** | Logarithmic or linear continuation coordinate. Units: logical. Values: .true.  or  .false.. | Scans. Log wavevector coordinates require positive endpoints; use linear for negative kpar. | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-scan_style"></a>`scan_style` (integer) | **No safe default** | Select family of scan coordinate. Units: enumeration. Values: -1 compound k;0 globals;1..nspec individual species. | Options1/2/4; legacy parameter scans. Mixed/cold allow wavevector scans only.  | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-ns"></a>`ns` (integer) | **No safe default** | Number of output intervals; normally ns+1 output points. Units: count. Values: positive for wavevector scans. | Scans. Internal steps ns\*nres; hidden substeps can fail without a written point. | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-nres"></a>`nres` (integer) | **No safe default** | Substeps per output interval. Units: count. Values: positive for wavevector scans. | Scans. Writes when mod(internal_index,nres)=0; companion point is internal index.  | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-heating"></a>`heating` (logical) | `.false.` | Request species absorption and channels; corrected also companion energy. Units: logical. Values: .true.  or  .false.. | Serial/2D, grid scan_input_2. Corrected collisionless wavevector scans only; second-axis flag in 2D. Maps reject corrected and warn/ignore legacy. Cold kinetic slots NaN, warm totals/channels retained. | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-heating_gamma_ratio_max"></a>`heating_gamma_ratio_max` (real) | `0.05` | Maximum abs(gamma/omega_r) for corrected heating applicability. Units: dimensionless fraction; 0.05 means 5%. Values: finite, nonnegative; zero permits only zero damping/growth under this check. | Resets to default for each &scan_input_N. Serial scans use their own setting; option=2 uses scan_input_2. Applicability bit 1 is set only when the ratio strictly exceeds this limit; equality passes. Relaxing it enables normalized Ps/channels only if all other checks pass; the independent energy-balance check is unchanged. Ignored by legacy heating. | Set `heating_gamma_ratio_max=0.05` in &scan_input_N. |
| <a id="param-eigen"></a>`eigen` (logical) | `.false.` | Request E/B and existing kinetic/cold density/velocity eigenfunctions. Units: logical. Values: .true.  or  .false.. | Serial/2D; grid uses scan_input_2. Serial uses this scan; 2D uses second axis. Map request warns and is ignored. Ex=1 field and per-species moment failures are separate.  | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |
| <a id="param-tensor"></a>`tensor` (logical) | `.false.` | Request six independent complex susceptibility elements per species. Units: ordinary chi. Values: .true.  or  .false.. | Serial/2D output; grid requires BOTH tensor flags. Serial scan flag; 2D requires both axes true. Map requests warn and are ignored. make_scan(..., tensor=False) exposes the independent option. | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in) |

<a id="fpc-fpc"></a>

## FPC: &fpc

Used by options 6/7 only, in all-kinetic configurations. Grids are historical species-parallel-thermal-speed coordinates; all component bounds and delv should be integer-compatible. FPC normalization and numerical moment closure have caveats; see the [FPC appendix](README-JETPLUME.md). 

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-vperpmin"></a>`vperpmin` (real) | **No safe default** | Velocity-grid minimum for vperp. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vperpmax"></a>`vperpmax` (real) | **No safe default** | Velocity-grid maximum for vperp. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input.  | See group workflow / tests |
| <a id="param-vparmin"></a>`vparmin` (real) | **No safe default** | Velocity-grid minimum for vpar. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vparmax"></a>`vparmax` (real) | **No safe default** | Velocity-grid maximum for vpar. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-delv"></a>`delv` (real) | **No safe default** | FPC velocity-grid spacing. Units: species parallel thermal-speed units in historical FPC grid. Values: positive. | FPC6/7. No comprehensive grid validation; zero spacing can cause invalid arithmetic. | See group workflow / tests |
| <a id="param-vxmin"></a>`vxmin` (real) | **No safe default** | Velocity-grid minimum for vx. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vxmax"></a>`vxmax` (real) | **No safe default** | Velocity-grid maximum for vx. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vymin"></a>`vymin` (real) | **No safe default** | Velocity-grid minimum for vy. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vymax"></a>`vymax` (real) | **No safe default** | Velocity-grid maximum for vy. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vzmin"></a>`vzmin` (real) | **No safe default** | Velocity-grid minimum for vz. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-vzmax"></a>`vzmax` (real) | **No safe default** | Velocity-grid maximum for vz. Units: species parallel thermal-speed units in historical FPC grid. Values: finite ordered velocity bounds. | FPC option6 gyro or7 Cartesian. All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input. | See group workflow / tests |
| <a id="param-elecdircontribution"></a>`elecdircontribution` (real) | `0.0` | Debug selection of electric forcing components for FPC. Units: selector. Values: 1 Ex only,2 Ey only,3 Ez only;0/other unchanged. | FPC6/7. Changes a debug diagnostic forcing, not dispersion. Not new mixed FPC support.  | See group workflow / tests |

<a id="experimental-radial-controls"></a>

## Experimental radial controls

The following &radial_input and selected &k_range fields are retained for option=5, not promoted to general mixed-plasma support. All kinetic only; corrected radial heating is rejected. No small standalone public radial model is currently validated as a primary example.

<a id="radial_input"></a>

### &radial_input

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-nrad"></a>`nRad` (integer) | **No safe default** | Number of radial model rows. Units: count. Values: positive required model-row count. | Experimental option5. Requires modelName0.rad and one species table per population.  | Experimental; no primary input |
| <a id="param-nrad_sub"></a>`nRad_sub` (integer) | `1` | Radial interpolation substeps. Units: count. Values: positive intended. | Option5. Default 1; interpolate in selected input variables, not generally invariant between normalizations. | Experimental; no primary input |
| <a id="param-modelname"></a>`modelName` (character(100)) | **No safe default** | Prefix of &lt;modelName&gt;0.rad globals and &lt;modelName&gt;s.rad species tables. Units: path. Values: valid input model basename/path. | Option5. Global rows radius,beta,vtp; species rows radius,tau,mu,alpha,Q,D,U.  | Experimental; no primary input |
| <a id="param-radial_heating"></a>`radial_heating` (logical) | **No safe default** | Request historical radial heating. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output). | Experimental; no primary input |
| <a id="param-radial_eigen"></a>`radial_eigen` (logical) | **No safe default** | Request radial fields/moments. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output).  | Experimental; no primary input |
| <a id="param-k_scan"></a>`k_scan` (integer) | **No safe default** | Radial wavevector selection:0 fixed,1 kz,2 kp,3 magnitude,4 angle,5 kp/kz plane,6 angle/magnitude plane. Units: enumeration. Values: 0..6. | Option5. Selects different &k_range namelist layouts; cold-containing radial unsupported.  | Experimental; no primary input |

<a id="k_range-choose-the-layout-matching-k_scan"></a>

### &k_range (choose the layout matching k_scan)

These controls belong only to the experimental radial pathway, not ordinary `scan_input_N`. The layout table below specifies which names may coexist. Angles are degrees; wavevector coordinates use the selected native length.

| Name / type | Default | Definition, units and allowed values | Applicability, dependencies and actual validation | Existing example |
|---|---|---|---|---|
| <a id="param-kperp_1"></a>`kperp_1` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout.| Experimental; no primary input |
| <a id="param-kpar_1"></a>`kpar_1` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout.  | Experimental; no primary input |
| <a id="param-kpar_2"></a>`kpar_2` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout. | Experimental; no primary input |
| <a id="param-nk"></a>`nK` (integer) | **No safe default** | Radial k_range output interval count for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate.  | Experimental; no primary input |
| <a id="param-kres"></a>`kres` (integer) | **No safe default** | Radial k_range resolution substeps for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate.  | Experimental; no primary input |
| <a id="param-rad_log"></a>`rad_log` (logical) | **No safe default** | Choose logarithmic interpolation for the named radial wavevector coordinate. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output).  | Experimental; no primary input |
| <a id="param-rad_heat"></a>`rad_heat` (logical) | **No safe default** | Request historical radial heating. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output).  | Experimental; no primary input |
| <a id="param-rad_eigen"></a>`rad_eigen` (logical) | **No safe default** | Request radial fields/moments. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output). | Experimental; no primary input |
| <a id="param-kperp_2"></a>`kperp_2` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout. | Experimental; no primary input |
| <a id="param-k_1"></a>`k_1` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout.  | Experimental; no primary input |
| <a id="param-k_2"></a>`k_2` (real) | **No safe default** | Initial/final wavevector coordinate in radial k_range layout. Units: native selected k-length convention. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout. | Experimental; no primary input |
| <a id="param-theta_1"></a>`theta_1` (real) | **No safe default** | Initial/final angle in radial k_range layout. Units: degrees. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout.  | Experimental; no primary input |
| <a id="param-theta_2"></a>`theta_2` (real) | **No safe default** | Initial/final angle in radial k_range layout. Units: degrees. Values: finite; positive if logarithmic. | Option5; k_scan selects active fields. Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout. | Experimental; no primary input |
| <a id="param-nk2"></a>`nk2` (integer) | **No safe default** | Radial k_range output interval count for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate.  | Experimental; no primary input |
| <a id="param-kres2"></a>`kres2` (integer) | **No safe default** | Radial k_range resolution substeps for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate. | Experimental; no primary input |
| <a id="param-rad_log_perp"></a>`rad_log_perp` (logical) | **No safe default** | Choose logarithmic interpolation for the named radial wavevector coordinate. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output).  | Experimental; no primary input |
| <a id="param-rad_log_par"></a>`rad_log_par` (logical) | **No safe default** | Choose logarithmic interpolation for the named radial wavevector coordinate. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output). | Experimental; no primary input |
| <a id="param-ntheta"></a>`ntheta` (integer) | **No safe default** | Radial k_range output interval count for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate. | Experimental; no primary input |
| <a id="param-thetares"></a>`thetares` (integer) | **No safe default** | Radial k_range resolution substeps for its designated coordinate. Units: count. Values: positive intended. | Option5, layout depends on k_scan. Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate. | Experimental; no primary input |
| <a id="param-rad_log_k"></a>`rad_log_k` (logical) | **No safe default** | Choose logarithmic interpolation for the named radial wavevector coordinate. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output). | Experimental; no primary input |
| <a id="param-rad_log_theta"></a>`rad_log_theta` (logical) | **No safe default** | Choose logarithmic interpolation for the named radial wavevector coordinate. Units: logical. Values: .true.  or  .false.. | Experimental radial option5. Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output). | Experimental; no primary input |

| k_scan | Required k_range layout |
|---|---|
| 0, fixed vector | kperp_1, kpar_1 |
| 1, fixed kperp / varying kpar | kperp_1, kpar_1, kpar_2, nK, kres, rad_log, rad_heat, rad_eigen |
| 2, fixed kpar / varying kperp | kperp_1, kperp_2, kpar_1, nK, kres, rad_log, rad_heat, rad_eigen |
| 3, fixed angle / varying magnitude | k_1, k_2, theta_1, nK, kres, rad_log, rad_heat, rad_eigen |
| 4, fixed magnitude / varying angle | k_1, theta_1, theta_2, nK, kres, rad_log, rad_heat, rad_eigen |
| 5, Cartesian plane | kperp_1, kperp_2, kpar_1, kpar_2, nk, kres, nk2, kres2, rad_log_perp, rad_log_par, rad_heat, rad_eigen |
| 6, polar plane | k_1, k_2, theta_1, theta_2, nk, kres, ntheta, thetares, rad_log_k, rad_log_theta, rad_heat, rad_eigen |

`nK` and `nk` are the same case-insensitive parameter. Model tables contain
radius in the first column, then beta/vtp in modelName0.rad or the six species
parameters in modelNameN.rad. nRad sets rows, nRad_sub interpolation substeps.
radial_eigen/radial_heating concern radial continuation; rad_eigen/rad_heat are
stored on the auxiliary wavevector-range scan. This historical two-level
interface has no primary curated example and should not be generalized without
separate validation. The existing synthetic radial test is the developer reference.

<a id="interactions-and-validation-gaps"></a>

## Interactions and validation gaps

- Cold-containing inputs require d_ref, collision_type=0, no FPC, and only
  supported wavevector scans. Cold references require betap=0; warm d_ref
  references require tauS_1=betap within the implemented relative comparison.
- Warm-reference adaptation needs positive reference beta, vtp and physical
  anisotropy. Cold anisotropy is ignored.
- Modern selectors and Bessel count/tolerance bounds are checked explicitly.
  Many old finite/count/storage/bound/Kn/FPC-grid checks remain incomplete.
- Self mass/charge/density ratios are a user consistency responsibility; the old
  mass-ratio warning is nonterminal and not a comprehensive reference check.
  Charge neutrality, zero current and reference drift are warning-only, not
  imposed equilibrium constraints.
- All-kinetic kperp=0 warns but remains unsupported. Every wavevector scan excludes
  kpar=0, including hidden points. Do not replace it with an arbitrary small floor. **This will be updated in a future release**
- Map eigen/tensor and legacy heating requests warn and remain ignored;
  corrected map heating is rejected. 
  See [Modes](modes.md) rather than assuming
  that a flag implies a diagnostic file exists.

<a id="python-wrapper-contract"></a>

## Python wrapper contract

`linfpclib.linfpc.plume_input` supports the scalar params/species/maps/scan/guess/
fpc groups above. It preserves low_n, both log flags, nr/ni, tensor and all
Cartesian FPC bounds on read/write roundtrips. `writeOut` is independent of
`use_map`, and `positive_roots` is independent of `loggridw`.

`set_maps(loggridw, omi, omf, gami, gamf, positive_roots,
loggridg=False, nr=128, ni=128)` retains its historical first six arguments.
`make_scan(..., tensor=False)` adds an optional tensor flag. Logical values
in wrapper dictionaries use canonical Fortran strings. Quoted paths, comments,
D exponents and scalar logical variants are supported.

`load_from_file` replaces prior indexed state, preserves an explicit dataName,
and rejects unknown fields/unsupported groups or array/expression syntax instead
of silently dropping controls. Assign `obj.dataname` **after loading** to redirect
output. `write_input` still takes an explicit outputname. Radial/k_range groups
are not supported by this wrapper. It is not a general Fortran namelist parser.
Examples and executable reader code are on [Output](output.md#python-workflows).
