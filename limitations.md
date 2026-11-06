<a id="limitations-and-troubleshooting"></a>

# Limitations and troubleshooting

[Manual home](./) · [Capability matrix](models.md#capability-matrix) · [Input](input.md)

This page consolidates current restrictions. A warning or accepted process exit
is not a statement that a calculation is physically supported. A valid dispersion
root also does not imply every supplementary diagnostic is available.

<a id="geometry-and-plasma-models"></a>

## Geometry and plasma models

- All-kinetic **exact parallel kperp=0 is unsupported**, in either normalization.
  Input/scan preparation warns; the historical kinetic evaluation may return
  invalid values. Do not use a tiny arbitrary kperp floor to claim the limit.
- Cold/mixed exact parallel transverse/kinetic limits are supported with nonzero
  kpar. Their successful tests do not certify an all-kinetic exact-axis dispatcher.
- Exact perpendicular kpar=0 is supported only for fully cold single-wavevector
  dispersion. Mixed/warm perpendicular response is deferred. Every wavevector
  scan, including fully cold scans, rejects kpar=0 at written or hidden points.
- Ex=1 cannot represent an Ex=0 mode. In particular the cold perpendicular
  ordinary root may exist but its electromagnetic eigenfunction is unavailable.
  No alternative field gauge is introduced. Poles, rank degeneracy and off-mode
  matrices can also invalidate fields.
- The equilibrium is nonrelativistic, homogeneous and gyrotropic, with parallel
  equilibrium drifts. Arbitrary distributions, perpendicular equilibrium drifts,
  relativity and MPI/OpenMP parallelization are not available.
- Cold-containing input requires collisionless d_ref. Cold anisotropy is ignored;
  there is no cold pressure or cold kinetic absorption/FPC. Positive tauS, however
  small, remains kinetic. Model-changing scans through zero are unsupported.

<a id="collisions-and-fpc"></a>

## Collisions and FPC

The all-kinetic collision_type=1 path is a nonconservative neutral-charged Krook
response. It is not a general Coulomb or energy-conserving collision operator.
The longitudinal collisional zz correction was independently validated and is
not undone by legacy numerical switches. A separate finite-rate xz/yz issue
remains outside the validated corrected-diagnostic scope. Do not generalize a
passing zz test to all collisional tensor elements.

Corrected heating requires collision_type=0. Explicit legacy all-kinetic
collisional diagnostics retain historical behavior, without a new accuracy claim.
Cold/mixed collisions are unsupported. Kn positivity/finiteness and some old
inputs are incompletely validated; a huge default Kn is not a general collision
model. Read [Input](input.md) and record actual collision_type and Kn.

Gyro/Cartesian FPC remains an all-kinetic historical pathway, including its thermal
solver/header conventions. There is **no mixed-plasma FPC support**, even for
warm species within a mixture. Susceptibility adaptive harmonics do not make
FPC's fixed internal harmonic reconstruction adaptive. Distribution moments are
debug integrations with substantial numerical limitations, not the validated
analytical species moments. Header phase-speed normalization has a reference-
anisotropy caveat; see [FPC appendix](README-JETPLUME.md).

Reader NaN preservation fixes interpretation, not FPC physics, normalization or
energy closure. Old species-index labels, header-offset retries and debug-only
readers still need scrutiny. Unavailable samples are not physical zeros.

The contour/resonance plotter's explicit **file reader** is not a universal
cold-reference overlay calculator. Its `vthpar_over_ref` helper still forms
`sqrt(mu*theta/beta_R)` in d_ref, which is undefined at beta_R=0. Reading a
cold-reference mode table can therefore work while a reference-thermal
resonance overlay does not. Use the physical units on Normalization and the
supported data readers; no cold-reference overlay extension is claimed here.

<a id="diagnostic-and-workflow-boundaries"></a>

## Diagnostic and workflow boundaries

- Corrected energy/heating is a real-frequency weak-mode approximation, not exact
  conserved energy for strongly growing/damped, degenerate or reactive modes.
  energy_status and applicability differ; formal H may be resolved while Ps is
  unavailable. Cold real-frequency poles can block energy for a valid complex mode.
- Map eigen/tensor requests and legacy map heating are ignored with warnings;
  corrected map heating is rejected. Map file schemas remain determinants and
  candidates. Use a supported root/diagnostic workflow to obtain moments or heat.
- Two-dimensional eigen/heating flags come from the second scan; tensor output
  needs both flags. Flag semantics are not identical across these diagnostics.
- Mixed/cold scans support wavevectors only. Historical plasma-parameter and
  radial workflows are all kinetic, with old restrictions and legacy heating.
  They are not automatically validated for arbitrary models or new diagnostics.
- Historical logarithmic **frequency maps** mix exponent origins with ratio
  increments and discard imaginary signs. They do not reach requested endpoints.
  Use linear map grids; [the exact equations](modes.md#linear-versus-historical-logarithmic-maps)
  remain unchanged pending a separate compatibility decision.
- nroot_max controls guesses/continued branch slots, not all .roots rows.
  Candidate files precede final filters. Sorting by least damping is not the
  enabled default. A mode label is a continuation slot, not proof of a universal
  physical branch identity through crossings.

<a id="diagnosing-a-run"></a>

## Diagnosing a run

| Symptom | Check and next action |
|---|---|
| Indexed namelist cannot be found | Use exact lowercase &species_N/&guess_N/&scan_input_N headers, alone on their lines; no trailing header comments; correct consecutive indices |
| Cannot open input/scratch | Use a .in suffix, writable input directory and a short path; the historical argument buffer is 100 characters |
| Cannot open output | Run from repo root or a deliberately prepared working directory; mkdir -p data/<dataName>; writeOut=false does not suppress opens |
| No .roots after supplied guesses | option1/nscan0 refines/prints but does not write map candidate files |
| No eigen/heating file after a map | Map diagnostics are ignored/rejected; a flag is not evidence of a supported writer |
| No .energy | Require corrected direct diagnostics or a supported scan with heating=true; no map/FPC/legacy energy companion |
| Too few/many roots | Inspect finite refined candidates, final filters and nroot_max branch slots; vary map bounds/grid or guesses; .roots is not a final capped list |
| NaN frequency at a scan point | Inspect root/status warning, poles and finite matrix; last valid seed is retained on guarded paths; reduce step or reinitialize physically |
| Finite frequency but NaN E/B | Ex=0 gauge, unresolved rank/nullvector, off-mode matrix or invalid response; check Lambda E and independent mode polarization |
| Only one species has NaN moments | Inspect its status, density and Doppler/cyclotron poles; D=0 gives undefined relative moments; valid other species/fields survive |
| Cold species has NaN Ps/channels | Expected unavailability, not zero heat; its H_s may still be valid |
| Warm Ps/channels NaN with finite H | Check applicability bits for strong/near-reactive mode, or species/normalizer status; do not force weak energy balance |
| H/fractions unavailable | wr=0, real pole, unresolved derivative, cancellation or absent field; no small H floor; retain valid root/moments |
| Bessel cap warning | Finite calculation continued with unresolved tails; increase nmax and compare tensors/roots, not just one small coefficient |
| Faddeeva nonrepresentable output | Lower-half-plane growth of Z can exceed binary64; distinguish range failure from a physical pole; no silent legacy fallback |
| Different rho/d_ref curves | Convert k using sqrt(beta_R alpha_R), vtp and all tau; compare identical physical scan paths, not just identical numbers |
| Missing values appear as zeros in a plot | Keep readers' default NaNs; visualization_nan_to_zero is an explicit plot-only opt-in |
| Reader gives plausible wrong diagnostics | Select S, eigen, heating, low_n and correct file family; parameters are not heating; do not use debug .eigen as .modeN |
| Ambiguous serial contour file | Call reader API with scan_index=N; CLI currently has no selector; 2D schema uses scan2 flags |
| Ragged low_n=false legacy file | Orphan Ew record is historical; use explicit handling or supported rectangular/corrected companion; never treat Ew as energy |
| Extra old mode files or mismatched companion counts | Use a new outputName/isolated directory; complete the run; companions and buffered modes are not atomic under interruption |

<a id="scientific-acceptance-is-more-than-process-success"></a>

## Scientific acceptance is more than process success

Check finite frequency, status, normalized determinant and matrix/nullvector
residuals. Compare species-current/continuity and Faraday identities when using
moments. Perform convergence checks in step size and harmonic cap/tolerance.
A very small raw determinant, a finite old Ew, a smooth continuation curve or a
plot without NaNs does not independently establish accuracy.

Validation gaps remain for some global count/storage/finite/Kn/warm-density
inputs and reference self-ratios. The wrapper is an intentionally restricted
scalar parser, not a general namelist editor. Its explicit rejection of unknown
controls prevents silent loss, but does not certify a physically meaningful run.
Preserve inputs, source revision and numerical options for reproducibility.
The later documentation/example/deployment stages will address remaining workflow
presentation issues; they do not silently extend these scientific capabilities.
