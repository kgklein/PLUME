# PLUME Input

This is a reference for the key input parameters used by PLUME.

## Plasma dispersion function backend

`&params` accepts `z_backend='faddeeva'` (default) or `z_backend='legacy'`.
Other values fail input validation. Old inputs without this setting select the
new accurate default; to reproduce historical collisionless results explicitly
select **both** `z_backend='legacy'` and `bessel_mode='fixed128'` (or
`bessel_mode='fixed', bessel_nfixed=128`). All Z/Bessel combinations are supported.
The previous collisional `zz` physics correction remains active in either mode.

The accurate backend evaluates $Z(\zeta)=i\sqrt{\pi}w(\zeta)$ with the Landau
analytic continuation. The public `zet_in(kpar,zin)` convention is unchanged:
positive `kpar` returns $Z(zin)$; negative `kpar` returns $-Z(-zin)$.
Kinetic-only, warm-reference mixed, auxiliary-normalized cold-reference mixed,
and weakly collisional calculations share this choice. Fully cold calculations
do not evaluate Z. This setting changes no normalization or harmonic cutoff.

The historical approximation has an accuracy defect near $|\zeta|=3$ that
cannot be repaired by tightening its iteration tolerance. Faddeeva removes that
defect but can cost more CPU time; compiled benchmarks and residual accuracy
limits are in the [backend report](tests/FADDEEVA_BACKEND_REPORT.md).
Bessel approximation, truncation and cancellation remain separate error sources.

For the accurate backend, nonfinite arguments, `kpar=0` calls, nonrepresentable
outputs and numerical failures produce assigned quiet NaNs. Dispersion guards
mark the evaluation invalid; root refinement rejects it and frequency maps mask
isolated invalid points. There are no floors, damping offsets, clamping or
fallback to legacy. A later valid evaluation resets the failure state.
In particular $\zeta=-14i$ is representable, whereas $\zeta=-30i$ exceeds
binary64 range. The standalone API distinguishes invalid input (1), zero `kpar`
(2), output range (3), and numerical failure (4), with 0 meaning success.
Production dispersion status retains existing conventions: 1 invalid frequency,
2 zero frequency, cold poles 2–4, and 5 nonfinite numerical response. The legacy
backend retains its historical failure behavior for reproducibility.

The Python wrapper's `set_params(..., z_backend='faddeeva')` reads/writes the
same setting; use `z_backend='legacy', bessel_mode='fixed128'` for compatibility.
See [accurate input](inputs/example/faddeeva_accuracy.in) and
[historical numerical input](inputs/example/legacy_numerics.in).

## Bessel harmonic controls

These optional `&params` settings apply to warm kinetic species in either
normalization, including warm-reference and auxiliary-normalized cold-reference
mixtures and the existing weakly collisional all-kinetic calculation.

| Input | Default | Meaning |
| --- | --- | --- |
| `bessel_mode` | `'adaptive'` | `'adaptive'` selects each species' cutoff; `'fixed'` uses `bessel_nfixed`; `'fixed128'` is an alias requiring `bessel_nfixed=128`. |
| `bessel_tol` | `1.e-10` | Positive ordinary and quadratic-weighted Bessel-tail bounds must both be below this value. Require finite `0 < bessel_tol < 1`. |
| `bessel_nmax` | `512` | Adaptive cutoff cap, independently applied to each warm species. |
| `bessel_nfixed` | `128` | Fixed maximum **absolute** harmonic order: sum `-N..N`, including zero. |
| `bessel_diagnostics` | `.false.` | Write per-species cache statistics to the `.bessel` sidecar. |

Both counts must be positive integers. The current 32-bit harmonic arithmetic
requires counts at most 46338 to avoid integer overflow in the existing `n*n`
expressions. Fixed mode is not limited by `bessel_nmax`. Adjacent coefficients
through `N+1` are stored separately; they do not add another susceptibility
harmonic. Exact parallel limits already supported by PLUME remain analytical.
No wavevector component or temperature is replaced by a floor.

The physical argument is
`lambda_s=(kperp*d_ref)^2*Qs^2*alphS*tauS_d/(2*muS)`.
In `rho_ref` it is `kperp_rho^2*Qs^2*alphS/(2*muS*tauS_rho*alphS_ref)`.
Thus different species generally require different cutoffs. The reference
anisotropy in the historical gyroradius is retained. Auxiliary normalization
changes neither this argument nor the physical input state.

Reaching `bessel_nmax` without satisfying both bounds **warns and continues**
at that cutoff. The warning identifies species, argument, tolerance, cap and
bounds. Identical repeated frequency evaluations do not repeat it. Metadata
retains unresolved-point history even if a subsequent scan point converges.
Nonfinite Bessel arithmetic or allocation/range failure is a numerical failure,
not a successful capped approximation.

`bessel_tol` bounds omitted Bessel weights, not the full susceptibility error:
resonances, cancellation, Bessel-coefficient accuracy, and (with legacy Z)
historical plasma-dispersion-function errors can dominate. Check each observable by tightening the tolerance and increasing
the cap; a substantially higher fixed count is another useful check. Fixed 128
is part of historical compatibility (also select `z_backend='legacy'`),
not an accuracy reference at large arguments. It keeps
the historical Bessel evaluator and compiler path; other counts use the new
vector Bessel calculation. Consequently differences from fixed 128 can also
include improvements in coefficient accuracy.

The Python wrapper's `set_params` accepts the same optional keyword arguments,
including `bessel_diagnostics=True`, and reads/writes them in input files.
See [the example](inputs/example/adaptive_bessel.in) and
[the technical report](tests/ADAPTIVE_BESSEL_CACHE_REPORT.md).

## Normalization conventions

`normalization='rho_ref'` is the default and preserves existing inputs.
`normalization='d_ref'` uses the inertial length of `species_1`. The reference
species must be first in either convention. The symbols $\mu_s=m_R/m_s$,
$Q_s=q_R/q_s$ (signed), $D_s=n_s/n_R$, and $\alpha_s=T_{\perp,s}/T_{\parallel,s}$
have the same meaning in both modes.

| Parameter | `rho_ref` | `d_ref` |
| --- | --- | --- |
| `kperp`, `kpar` | $k\rho_R$, $\rho_R=v_{t\perp,R}/|\Omega_R|$ | $kd_R$, $d_R=c/\omega_{pR}$ |
| `vtp` | $v_{t\parallel,R}/c$ | $v_{A,R}/c$ |
| `betap` | $8\pi n_RT_{\parallel,R}/B_0^2$ | same, and equals reference `tauS` |
| `tauS` | $T_{\parallel,R}/T_{\parallel,s}$ | $2T_{\parallel,s}/(m_Rv_{A,R}^2)$ |
| `muS`, `Qs`, `Ds`, `alphS`, `vvS` | $m_R/m_s$, $q_R/q_s$, $n_s/n_R$, $T_{\perp,s}/T_{\parallel,s}$, $U_s/v_{A,R}$ | same |
| `Kn` | $(n_0\sigma\rho_R)^{-1}$ | $(n_0\sigma d_R)^{-1}$ |

For $\alpha_R=T_{\perp,R}/T_{\parallel,R}$ and positive reference temperature,
$k d_R=k\rho_R/\sqrt{\beta_{\parallel,R}\alpha_R}$,
$v_{t\parallel,R}/c=\sqrt{\beta_{\parallel,R}}\,v_{A,R}/c$, and
$\tau_s^{(d)}=\beta_{\parallel,R}/\tau_s^{(\rho)}$. In `d_ref`,
$v_{t\parallel,s}^2/v_{A,R}^2=\mu_s\tau_s^{(d)}$,
$\beta_{\parallel,s}=D_s\tau_s^{(d)}$,
$\omega_{ps}^2/\Omega_R^2=D_s\mu_s/[Q_s^2(v_{A,R}/c)^2]$,
and $\Omega_s/\Omega_R=\mu_s/Q_s$. The signed charge ratio is retained in
the cyclotron frequency and Doppler shift.

`d_ref` requires `species_1` `tauS=betap`.
A scan of either reference beta or reference `tauS` updates both.
Its other inputs remain fixed during that scan;
this is a different physical scan from changing only `betap` in `rho_ref`.
Two-dimensional scans cannot vary both redundant quantities independently.
Positive `tauS` selects a kinetic species. In collisionless `d_ref` runs,
exactly zero `tauS` selects a cold-fluid species; negative temperatures
are invalid. A cold reference requires `betap=0`; additional species may be warm or cold.
See the matched example files in `inputs/example/normalization_*.in`.

### Mixed kinetic/cold calculations (Stage 2)

Use `normalization='d_ref'`, `collision_type=0`, a warm first species, and
`nscan=0` for single-wavevector runs (wavevector scans are described below). Additional species may have exactly `tauS=0`; every positive
temperature, however small, remains kinetic. Cold `alphS` is ignored (zero is
allowed). Density, signed charge, mass and parallel drift retain their definitions.
No cold temperature is replaced by a numerical floor.
Parallel (`kperp=0`) and oblique propagation are supported. Mixed runs require
nonzero `kpar`: the warm susceptibility's exact perpendicular limit is not
implemented. No small nonzero wavevector is substituted automatically.

- `option=-1` evaluates the determinant and reports its status. For mixed runs
  Stage 4B.1/4B.2 also prints E/B and species moments at a validated
  eigenfrequency. Corrected collisionless heating/wave energy is available in
  Stage 4B.3A; FPC remains unavailable for cold-containing plasmas.
- `option=1`, `use_map=.false.` refines supplied `&guess_N` values with the
  existing secant method. Set `writeOut=.true.` to print the solutions.
- `option=0` searches a frequency map. Singular samples remain NaN in map
  output and are excluded from root selection; failed refinements are skipped.
- Cold-containing plasma-parameter/radial scans and FPC
  are unavailable; corrected collisionless heating is available; collisionless wavevector scans are supported in Stage 4A. Requests fail clearly; temperature scans that reach or leave
  zero are also rejected. All-kinetic scan and diagnostic behavior is retained.

Cold species with `rho_ref` and cold-containing collisional configurations
are unsupported. The kernel's zero-frequency,
Doppler and cyclotron poles are not regularized. A supplied guess that reaches
a pole or fails to converge reports an error rather than a fabricated solution.

Charge neutrality and zero equilibrium current are **warning-only** checks,
emitted during input reading. Densities/drifts are never adjusted. Nonstandard
homogeneous equilibria may require external forces or currents to maintain them.

See `inputs/example/mixed_warm_protons_cold_electrons.in` and the
[Stage 2 implementation report](tests/MIXED_DISPERSION_REPORT.md).

### Fully cold calculations (Stage 3A)

Set `normalization='d_ref'`, `collision_type=0`, `betap=0`, and `tauS=0`
for **every** species, including species 1. Supply positive `vtp=v_A,ref/c`.
The inertial length `d_ref=c/omega_p,ref` remains finite at zero temperature;
no thermal speed or gyroradius is constructed. Thus `rho_ref` cannot be used
for this configuration. Cold `alphS` has no physical meaning and is ignored:
zero, unity, or other finite placeholders give identical results.

| Configuration | Supported |
| --- | --- |
| `d_ref`, all species cold, `betap=0` | Yes, Stage 3A |
| `d_ref`, warm reference, other warm/cold species | Yes, Stage 2 |
| `d_ref`, cold reference, warm/cold additional species, `betap=0` | Yes, Stage 3B.2 |
| `rho_ref`, all species warm | Yes, unchanged legacy calculation |
| `rho_ref`, any cold species | No |

Fully cold runs accept exact `kperp=0`, exact `kpar=0`, oblique wavevectors,
and arbitrary signed parallel drifts `vvS`. The mixed-run restriction on exact
perpendicular propagation does not apply. Frequencies remain `omega/Omega_ref`;
positive imaginary frequency means growth. No small wavevector, temperature,
frequency offset, or artificial damping is inserted.

Use the same `option=-1`, `0`, or `1` root/determinant workflows described above,
with `nscan=0` for a single wavevector, or Stage 4A wavevector scans below. The existing secant and map search are reused. Root guesses
select branches and must avoid exact poles; map resolution must resolve the
branch being sought. Corrected collisionless heating/wave energy is available
(see Stage 4B.3A below); FPC, cold collisions, and plasma-parameter/radial scans remain unavailable.
Warning-only charge/current checks and their external-maintenance caveat apply
also to fully cold equilibria.

Run `./plume.e inputs/example/cold_proton_electron.in` from the repository root.
This example selects the parallel L (Alfvén/ion-cyclotron) branch at
`kpar*d_ref=0.06`. The companion
`inputs/example/cold_proton_electron_perpendicular.in` selects the low-frequency
X (fast magnetosonic) branch at exactly `kpar=0`, `kperp*d_ref=0.2`.
Both retain finite electron mass and print their roots to stdout.
The Python wrapper accepts the same inputs through
`set_params(..., normalization='d_ref')` and `make_species(tauS=0, ...)`.
See the [Stage 3A implementation report](tests/FULLY_COLD_DISPERSION_REPORT.md)
for tensor, root, singularity and regression validation.

### Cold-reference mixtures with warm species (Stage 3B.2)

Use `normalization='d_ref'`, `betap=0`, `tauS=0` for species 1,
`collision_type=0`, and `nscan=0` for single-wavevector runs (wavevector scans are described below). Additional species use exactly zero `tauS`
for cold-fluid response or positive `tauS` for kinetic response; negative values
are invalid. No temperature threshold is used. Supply positive `vtp=v_A,ref/c`.
The physical reference mass, signed charge, density, gyrofrequency, Alfvén
speed and inertial length stay fixed, as do species order and equilibrium drifts.
Cold anisotropy placeholders are ignored; warm anisotropies remain physical.

For oblique warm responses, PLUME selects the first warm species `j` and sets
`b_star=tauS(j)`. This internal numerical policy adds no required input parameter.
It constructs `V=vtp*sqrt(b_star)`, `K=k*d_ref*sqrt(b_star)` and inverse ratios
`t_s=b_star/tauS(s)`, with a separate normalization anisotropy of unity.
The validated kinetic expressions then give the same physical susceptibility,
independent of this arbitrary scale up to finite-precision conditioning.
Neither the physical reference temperature nor any species temperature is
replaced: `betap` and reference `tauS` remain exactly zero. No fictitious
population, reference thermal speed or gyroradius is introduced.

Exact parallel (`kperp=0`) propagation uses the existing analytical kinetic
limit with `v_t,s/v_A=sqrt(muS*tauS)`. Oblique propagation, either sign of `kpar`,
parallel drifts of warm/cold species, and complex frequencies in either
half-plane are supported. **Exact `kpar=0` with any warm species is unsupported.**
Use `option=-1` for direct evaluation, `option=1`, `use_map=.false.` for supplied
guesses, or `option=0`/`use_map=.true.` for supported frequency-map searches.
Set `writeOut=.true.` to print roots. The existing secant method is unchanged;
nearby branch-specific guesses are useful when several physical modes coexist.

Extreme temperature spreads or auxiliary scales can produce unrepresentable
intermediate kinetic products even though the normalization is algebraically
invariant. Invalid or nonfinite evaluations return the existing failure status
and NaN tensors, with no added damping or floors. A subsequent valid trial is
recomputed. The first-warm-species policy is provisional, not an optimized
conditioning guarantee. Species reordering can change that scale but must not
change the physical result beyond the documented numerical tolerance.

Cold-containing FPC, collisions,
plasma-parameter/radial scans and model-transition scans remain unavailable.
Wavevector scans are supported as described below. Charge
neutrality and current balance remain warning-only checks made at input reading;
nonstandard homogeneous equilibria may require external currents or forces.

Run `./plume.e inputs/example/mixed_cold_protons_warm_electrons.in` from the
repository root for an oblique damped mode with cold protons and warm electrons.
The Python wrapper accepts the same convention and restrictions. See the
[Stage 3B.2 implementation report](tests/COLD_REFERENCE_MIXED_REPORT.md) for
interfaces, cache behavior, independent benchmarks and numerical limits.

## Namelists in input files.

The following namelists and associated input parameters are read in by PLUME from the input file.

### *&params*  
General system parameters.

**`normalization`**
`'rho_ref'` (default) or `'d_ref'`.
Invalid values stop with a diagnostic.

**`betap`**  
Initial reference parallel plasma beta $\beta_{ref,\parallel} = 8 \pi n_{ref} T_{ref,\parallel}/B^2$.

**`kperp`**  
Initial perpendicular wavevector $k_{\perp} \rho_{ref}$ or $k_{\perp}d_{ref}$,
according to `normalization`.

**`kpar`**  
Initial parallel wavevector $k_{\parallel} \rho_{ref}$ or $k_{\parallel}d_{ref}$.

**`vtp`**  
Initial reference parallel thermal velocity $v_{t,ref,\parallel}/c$ in `rho_ref`,
or reference Alfvén speed $v_{A,ref}/c$ in `d_ref`.

**`nspec`**   
Number of plasma species or components.

**`nscan`**   
Number parameter or wavevector scans to execute.

**`collision_type`**   
Type of collision for the system.
- 0: no collisions; standard PLUME
- 1: non-conservative Krook Operator
     neutral-charged collisions only.

**`Kn`**   
Neutral Knudson number over characteristic scale $(n_0 \sigma \rho_{ref})^{-1}$.
or $(n_0 \sigma d_{ref})^{-1}$.
Ignored for `collision_type`=0 (collisionless evalulation).

**`option`**   
Determines set of scans to perform. Choice of:

- -1: Calculate the value of the dispersion relation
      at a single $(\omega/\Omega_{ref}, \gamma/\Omega_{ref})$
-  0: Calculate roots for input plasma parameters.
-  1: Calculate roots for input plasma parameters or Reads in root value
      and then scan over plasma parameters, with range and type specified in *.in file.
-  2: Calculate roots for input plasma parameters or Reads in root value
      and then scan over two-dimensional plasma parameters space,
      with range and type specified in *.in file.
-  3: Deprecated...
-  4: Make multiple maps of complex frequency space.
-  5: Find roots for parameters along a prescribed path
      Path is set by solar wind models, with values calculated and
      output by helper function (in development, the radial scan function.).

**`nroot_max`**   
Number of dispersion solutions to find and follow.

**`use_map`**   
Choice of:  

- True: Searching for roots over a map in complex frequency space (see &maps namelist).  
- False: Input `nroot_max` guesses for solutions (see &guess_1 namelist).

**`low_n`**  
Logical to toggle on or off outputing the $n=0$ and $\pm 1$ resonances for
characterizing heating channels. Separates Landau, Transit, and Cyclotron damping.
This follows [Huang et al 2024 JPP](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174) and is expressed in Sec. 5 of [Klein, K. G., Howes, G. G.,
and Brown, C. R., 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2).

**`new_low_n`**  
*Depricated* Logical to toggle on or off outputing the $n=0$ and $\pm 1$ resonances for
characterizing heating channels. Previous had multiple versions of the low n susceptability calculation. As we are only using Eqns. 39,40, and 41 from [Klein, K. G., Howes, G. G.,
and Brown, C. R., 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2), this logical flag is unnecessary.

**`writeOut`**  
Verbosity control for screen output, independent of `use_map`. It does not
disable file output or required output directories. Input warnings remain
visible with `writeOut=.false.`. The Python wrapper preserves both switches
independently.

**`dataName`**  
Subdirectory (below PLUME/data/) where outputs will be written.

**`outputName`**
String to be included in output files to identify specific run.

### *&guess_m*  
Initial guess of complex frequency for $m$th solution.  
Only used when `use_map`=.false.  
Need to have number of namelists equal to `nroot_max`.

**`g_om`**  
Guess for real solution $\omega_{r}/\Omega_{ref} $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_{ref} $.


### *&maps*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{r}/\Omega_{p}$ map search.
Only linear grids use `omi` and `omf` as physical endpoints; see the
historical logarithmic convention below.

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{p}$ map search.
Only linear grids use `gami` and `gamf` as physical endpoints.

#### Historical logarithmic maps

The logarithmic map equations are retained for compatibility, but their
origins and increments use inconsistent conventions. For i=0..nr, j=0..ni:

```
wr(i) = 10**(omi + i*(log10(omf)-log10(omi))/nr)
wi(j) = -10**(abs(gami) + j*(log10(abs(gamf))-log10(abs(gami)))/ni)
```

Real logarithmic bounds must be positive; imaginary magnitudes must be
nonzero. Both imaginary input signs are discarded and all logarithmic
imaginary nodes are negative. Thus `omi=.1, omf=1` gives approximately
1.258925 to 12.58925; `gami=-1e-4, gamf=-1e-2` gives approximately
-1.000230 to -100.0230. These are neither physical endpoint bounds nor a
consistent exponent-bound convention. PLUME and the Python wrapper warn
when either log flag is requested. Use linear grids for physical endpoint
control. Stage 5B.0 does not change these equations; a future correction
requires an explicit compatibility decision. See the
[characterization report](tests/PLUME_5B0_CORRECTIONS.md).

**`omi`**  
Initial $\omega_{r}/\Omega_{p}$ endpoint for a linear complex map search;
logarithmic behavior is defined above.

**`omf`**  
Final $\omega_{r}/\Omega_{p}$ endpoint for a linear complex map search.

**`gami`**      
Initial $\gamma/\Omega_{p}$ endpoint for a linear complex map search.

**`gamf`**  
Final $\gamma/\Omega_{p}$ endpoint for a linear complex map search.

**`nr`**  
Number of $\omega_{r}/\Omega_{p}$ intervals (`nr+1` nodes) in the frequency grid.
Default is 128 if not specified.

**`ni`**  
Number of $\gamma/\Omega_{p}$ intervals (`ni+1` nodes) in the frequency grid.
Default is 128 if not specified.

**`positive_roots`**
If true, the map routine will only select solutions with $\omega_{|textrm{r}} \ge 0$.

### *&species_j*  
Species or component parameters list for distribution $f_{j}$.

**`tauS`**  
Relative parallel temperature $T_{ref,\parallel}/T_{j,\parallel}$ in `rho_ref`;
Alfvénic temperature $2T_{j,\parallel}/(m_{ref}v_{A,ref}^2)$ in `d_ref`.
Exactly zero selects cold-fluid susceptibility under the Stage 2/3A/3B.2 restrictions
above; it is not a valid legacy inverse-temperature ratio.

**`muS`**  
Relative mass $m_{ref}/m_{j}$.

**`alphS`**  
Temperature anisotropy $T_{j,\perp}/T_{j,\parallel}$; ignored for cold species.

**`Qs`**  
Relative charge $q_{ref}/q_{j}$.

**`Ds`**  
Relative density $n_{j}/n_{ref}$.

**`vvS`**  
Drift speed parallel to mean magnetic field, normalized to the reference Alfven velocity
$v_{j,drift}/v_{A,ref}$, where $v_{A,ref}=B/\sqrt{4 \pi n_{ref} m_{ref}}$

**`nunS`**  
Neutral-charged species collision frequency, $\nu_{n,s}/\Omega_{ref} = n_0$


## Mixed wavevector scans (Stage 4A)

All supported compositions use the existing scan framework and `disp/rtsec`:
all kinetic (`rho_ref` or `d_ref`), fully cold (`d_ref`), warm-reference mixed
(`d_ref`), and cold-reference mixed (`d_ref`, `betap=0`). Any number of supported
warm and cold species, parallel drifts and signed nonzero `kpar` are permitted.
Cold-containing scans must be collisionless. No species can change model.

| Workflow | Inputs | Traversal |
| --- | --- | --- |
| Serial continuation | `option=1`, `nscan>=1` | Each scan starts from the preceding scan's final valid root. |
| Two-dimensional grid | `option=2`, `nscan=2` | Outer-axis continuation; each inner traversal starts from its saved row-start seed. |
| Repeated frequency maps | `option=4`, `nscan=1` | Independent existing map search at each point; no branch association between maps. |

Wavevector coordinates use the existing namelists:

| `scan_style` | `scan_type` | Existing interpretation |
| ---: | ---: | --- |
| 0 | 0 | `kperp`, from `swi` to `swf` |
| 0 | 1 | `kpar`, from `swi` to `swf` |
| -1 | 0 | Current vector to `(kperp,kpar)=(swi,swf)` |
| -1 | 1 | Current `atan(kperp/kpar)` to `swi` degrees, fixed signed magnitude; `swf` unused |
| -1 | 2 | Current signed magnitude to `swf`, fixed angle; `swi` unused |

The historical signed magnitude is `kpar/cos(atan(kperp/kpar))`; for negative
`kpar`, use a negative target `swf` to retain that propagation direction.
Logarithmic scans require positive scanned coordinates; use linear scans for
negative parallel coordinates/signed magnitudes. Exact parallel `kperp=0`
uses the existing mixed/cold dispersion limits. Geometry support never extends
a susceptibility limit that the single-point evaluator does not implement.

Two-dimensional wavevector grids support distinct Cartesian `kperp/kpar` axes,
or ordered angle (`scan_input_1`) and magnitude (`scan_input_2`) axes. Overlapping
coordinate definitions are rejected. Existing all-kinetic plasma-parameter
scans are retained, but no plasma-parameter or radial scans are enabled for
cold-containing configurations in Stage 4A.

Require positive `ns` and `nres`; the existing scan takes `ns*nres` continuation
steps and writes `ns+1` rows. Input bounds must be finite. Scans are preflighted
before root initialization. If **any** sampled point has `kpar=0`, including a
90-degree angle endpoint, the run stops with a Stage 4A diagnostic. No epsilon
wavevector is inserted. This restriction includes fully cold scans but does
not change fully cold single-wavevector perpendicular calculations.

Cold-containing runs initialize their roots at the actual first scan point.
Their physical plasma parameters and reference definition are fixed throughout:
wavevectors are `k*d_ref`, frequencies are `omega/Omega_ref`, and `vtp=v_A/c`.
Warm-reference `d_ref` conversions retain `sqrt(betap*alphS_ref)`; cold-reference
runs bypass the thermal adapter and retain `b_star=tauS(first warm species)`.
Changing only `kpar` reuses Bessel coefficients; changing `kperp` updates/resizes
the existing cache using the physical argument. Z/Bessel backend choices remain
independent and unchanged.

For cold-containing scans `eigen=.true.` enables
E/B fields and species density/velocity moments. Failed species moment blocks
are marked NaN. Both flags default false.
Heating requests require `heating_convention='corrected'` and collisionless inputs.
Optional `tensor=.true.` uses the existing susceptibility output and introduces
no new columns. Stage 4B.2 supplies density/velocity perturbations;
FPC for mixed plasmas remains outside
this stage, including the agreed later policy for kinetic species. Cold values
are unavailable, not physical zeros. All-kinetic diagnostics remain available.

A failed root is written as NaN in the existing frequency columns; tensor and
requested diagnostic slots are invalidated as well. A warning gives
the wavevector and failure flag. Only successful roots update the continuation
seed. Later points retry from the last valid seed. An intermediate point hidden
by `nres` can fail and warn even if no NaN output row is scheduled there.
No branch tracker is introduced: reduce steps and inspect crossings, poles,
strong damping and recovery after failures. Cold-containing map scans retain
the supplied `&maps` frequency bounds, avoiding undefined cold anisotropy;
legacy all-kinetic map-span behavior is unchanged.

Examples: [warm reference](inputs/example/mixed_scan_warm_reference.in),
[cold reference](inputs/example/mixed_scan_cold_reference.in), and
[fully cold](inputs/example/fully_cold_scan.in). The Python wrapper accepts the
same wavevector scans through `make_scan`; its cold-containing writer validation
rejects parameter scans, mixed FPC, and legacy mixed-heating requests.

## Electromagnetic eigenfunctions (Stage 4B.1)

Supported fields use the complete production `lam` matrix, retaining each
species contribution in `susc(species,row,column)`. Existing kinetic susceptibility,
cold susceptibility, auxiliary normalization and numerical backend selections
are unchanged. Eigenfunctions are right nullvectors of the ordinary complex
wave matrix, not Hermitian eigenvectors at complex frequency.

For cold/mixed scans use `eigen=.true.`; corrected collisionless heating may also be enabled. Serial and 2D Stage 4A
workflows and both signs of nonzero kpar are supported. Option 2 retains its
existing rule that **scan_input_2 selects the output diagnostics**. Stage 4B.2
populates the existing species velocity/density slots using cold-fluid or kinetic
moments as appropriate. All-kinetic calculations retain their full existing
moments behavior. Corrected kinetic heating in mixtures is available in Stage
4B.3A; FPC remains deferred to Stage 4B.3B.

`option=-1` prints fields and species moments for a cold-containing input when the supplied
frequency has a validated electromagnetic nullvector. It prints an unavailable
status for off-mode or invalid frequencies. Root-only `option=1, nscan=0`
retains its existing root-only output; use a known root with option=-1, an
existing fixed-endpoint eigen scan, or the Fortran API for fields.

The amplitude/phase gauge is Ex=1+0i. The fields are `(Ex,Ey,Ez)` and `(Bx,By,Bz)`
with B0 along z and k=(kperp,0,kpar). The Fourier convention is exp(i*k.r-i*omega*t).
No handedness or component ordering changes. In Gaussian field units both E
and B outputs are divided by the physical deltaEx amplitude. To convert B to
units of B0, multiply by deltaEx/B0. These are not unit-energy eigenfunctions.

Faraday in PLUME's normalized inputs is

- `omega_hat * B = cross(k_rho/(vtp*sqrt(alphS_ref)), E)` for rho_ref;
- `omega_hat * B = cross(k_d/vtp, E)` for d_ref, where vtp=v_A/c.

The frequency is full complex omega/Omega_ref, including damping/growth, and
not just its real part. The anisotropic rho_ref conversion gives
`k_rho/(vtp_solver*sqrt(alphS_ref))=k_d/(v_A/c)`. No beta or auxiliary b_star
factor remains. Cold-reference fields never require a reference thermal speed.

The existing two-row cofactor elimination is reused. Equivalent row-pair
cofactors handle zero fixed divisors, including exact parallel mixed/cold
propagation. No LAPACK/general eigensolver or regularization is introduced.
Numerically rank-deficient matrices and modes with vanishing/unresolvable Ex
are unavailable in this gauge. The exactly perpendicular cold O mode is one
such case; the perpendicular X mode has a valid Ex gauge and is supported for
single-wavevector calls. Stage 4A still rejects perpendicular scan points.
Existing all-kinetic propagation limits are unchanged.

Field-only Fortran interfaces in `electromagnetic_fields`:

```fortran
! Refresh the total production matrix through the existing callback:
call calc_em_eigen(omega, electric, magnetic, ierr, disp)
! Or calculate from an explicit complete matrix, without global state:
call calc_em_fields(matrix, omega, kpar, kperp, vtp, normalization, &
                    alpha_ref, electric, magnetic, ierr)
```

Status 0 means success, 1 invalid inputs/evaluation, 2 degenerate/ill-conditioned
rank, 3 unavailable Ex gauge, 4 no acceptable nullvector at this frequency,
and 5 nonfinite arithmetic. Failed fields contain assigned quiet NaNs.
Field failure does not turn a valid dispersion matrix into an invalid root.
The direct matrix API accepts native input units and ignores alpha_ref for d_ref.
The historical `disprels::calc_eigen` remains the all-kinetic moment/heating API.

The historical `Ew` slots do not export the normalizer. Stage 4B.3A provides
complete corrected energy in a separate `.energy` companion; neither convention
is an exact conserved energy for arbitrary growing/damped modes. See below.

## Species moments (Stage 4B.2)

For `eigen=.true., heating=.false.`, cold-containing single-point and serial/2D
wavevector workflows now calculate cold-fluid moments for exactly `tauS=0`,
and kinetic moments for `tauS>0`. Negative temperature remains invalid. Species
ordering and physical density/mass/signed charge/drift definitions are unchanged.
No new input switch is required. Both signs of nonzero kpar, parallel/oblique
propagation, drifts and full complex frequencies are supported under the existing
dispersion restrictions. Exact perpendicular scans remain outside Stage 4A.

With `epsilon=deltaEx/B0`, the velocity columns are
`v=deltaU/(c*epsilon)` and density is `eta=(delta n/n0)/epsilon`.
They are **not** deltaU/vA or the absolute density delta n/n_ref.
In d_ref define `a=vtp=vA/c`, `kappa=k*d_ref`, `u=vvS`, and
`delta=omega_hat-kappa_parallel*u`. The cold momentum and continuity equations are

```text
-i*delta*v = (muS/Qs)*(E + v cross zhat + a*u*zhat cross B)
delta*eta = (kappa dot v)/a
```

The warm calculation reuses the existing species current from chi, removes the
`a*u*eta*zhat` equilibrium-drift density current, and uses continuity. It needs
no thermal reference conversion after chi is evaluated. For rho_ref equivalence
use `a=vtp/sqrt(betap)` and `kappa=k_rho/sqrt(betap*alphS_ref)`. The legacy
all-kinetic output/arithmetic remains unchanged. Neither cold nor warm physical
moments depend on the cold-reference auxiliary normalization.

The explicit-state module `species_moments` provides:

```fortran
call calc_cold_moments(species, omega, kpar_d, kperp_d, va_ref_c, &
                       E, B, velocity, density, status)
call calc_kinetic_moments(species, omega, kpar_d, kperp_d, va_ref_c, &
                          E, chi, velocity, density, status)
```

Both return one species in existing output units. Cold moments are derived from
momentum/continuity and do not call susceptibility, Bessel or Z routines, or read
temperature/anisotropy. The mixed adapter offers optional per-species status:

```fortran
call calc_mixed_eigen(omega,E,B,U,n,P,Psplit,.false.,field_status,disp, &
                      moment_status)
```

`field_status` retains its Stage 4B.1 meaning. Each `moment_status` is 0 valid,
1 invalid input, 2 zero frequency, 3 exact cold Doppler pole, 4 exact cold
cyclotron pole, 5 unrepresentable/nonfinite arithmetic, or 6 absent zero-density
species. A zero-density species has no defined relative-density/bulk-velocity
moment. Failure sets only that species' moment block to quiet NaNs and emits a
species/status warning, leaving valid fields and other species intact. Failed
fields leave all moments unavailable. A subsequent valid point is recomputed.
No floors, damping or offsets regularize poles. The warm current elimination
avoids a removable intermediate division at delta=0 when chi itself is finite;
it does not remove any cold response pole.

The existing column order, phase, scientific precision, and root continuation
are preserved. Metadata identifies model/units/status meanings. This moment stage does not change FPC, pressure moments or propagation.
Corrected heating/wave energy was added separately in Stage 4B.3A below.
See [validation and derivation](tests/COLD_FLUID_EIGENFUNCTION_REPORT.md).

### *&scan_input_l*
Inputs for scanning parameter space for $l$th scan.  

**`scan_type`**  
Type of parameter scan.

For `scan_style`=-1 (Global Two-Component Scan), options of:

- 0: Scan from $\textbf{k}_0 \rho_{ref}$ to $\textbf{k}_1 \rho_{ref}$.
Scans from current value of $(k_\perp,k_\parallel) \rho_{ref}$ to $k_\perp \rho_p$=`swi` and $k_\parallel \rho_{ref}$=`swf`.
- 1: Scan from $\theta_0$ to $\theta_1$ with fixed $|k|\rho_{ref}$.
Scans from current value of $(|k|\rho_{ref},\theta)$ to $(|k|\rho_{ref},$`swi`$)$, with `swi` in degrees.
- 2: Scan from $|k|_0 \rho_{ref}$ to $|k|_1 \rho_{ref}$ with a fixed value of $\theta$.
Scan from current value of $(k_\perp,k_\parallel) \rho_{ref}$ to $|k|\rho_{ref}$=`swf` with constant $\theta = \atan (k_\perp/k_\parallel)$.

*NOTE:* If you are executing a two-parameter plane scan of $\theta$ and $|k|\rho_{ref}$, $\theta$ (`scan_style`=-1, `scan_type`=1) must be set to `scan_input_1` and $|k| \rho_{ref}|$ (`scan_style`=-1, `scan_type`=2) to `scan_input_2`.
The opposite order will not properly scan through $|k| \rho_{ref}|$.
See
```
inputs/example/example_map_perp.in
```
for an example.

For `scan_style`=0 (Global Scan), options of:

- 0: $k_\perp \rho_{ref}$
- 1: $k_\parallel \rho_{ref}$
- 2: $\beta_{ref,\parallel}$
- 3: $v_{t,ref,\parallel}/c$

`swi` and `swf` represent the initial and final scan values for `scan_style`=0.

For `scan_style`=1-`nspec` (Parameter scan of component `scan_style`), options of:

- 0: $T_{ref,\parallel}/T_{j,\parallel}$.
- 1: $m_{ref}/m_{j}$.
- 2: $T_{j,\perp}/T_{j,\parallel}$.
- 3: $q_{ref}/q_{j}$.
- 4: $n_{j}/n_{ref}$.
- 5: $v_{j,drift}/v_{A,ref}$.

`swi` and `swf` represent the initial and final scan values for `scan_style`=1 through `nspec`.

**`scan_style`**  
Class of parameter scan. Options of:

- -1: Global Two-Component Scan
- 0: Global Scan
- 1-nspec: Single-Parameter Scan

**`swi`**  
Scan variable to define start of scan through parameter space.

**`swf`**  
Scan variable to define end of scan through parameter space.

**`swlog`**  
Use $\log_{10}$ (T) or linear (F) spacing.

**`ns`**  
Number of output scan values.

**`nres`**  
Resolution between output scan values.

**`heating`**  
Calculates heating rates if true.

**`eigen`**  
Calculates eigenfunctions if true. In cold-containing runs, Stage 4B.1 supplies
E/B fields and cold-fluid/kinetic species density and velocity moments;
failed moment blocks are unavailable (NaN).

**`tensor`**  
Outputs susceptibility tensor if true.     


## Heating and wave energy (Stage 4B.3A)

In `&params`, `heating_convention='corrected'` (default) or `'legacy'` selects
only the diagnostic convention. It does not select a Z backend, Bessel mode,
or normalization. For the new convention use `collision_type=0`.

```fortran
&params
heating_convention='corrected'
z_backend='faddeeva'
bessel_mode='adaptive'
! other ordinary plasma/root parameters ...
/
&scan_input_1
! existing wavevector scan coordinates ...
eigen=.true.
heating=.true.
/
```

`option=-1` computes diagnostics at the supplied frequency only when it has
an acceptable Ex=1 electromagnetic eigenfunction. With `option=1`/`option=2`,
heating follows the existing scan switches (in 2D, `scan_input_2` controls
outputs). Fully kinetic, fully cold, warm-reference mixed and cold-reference
mixed collisionless plasmas are supported. Existing geometry restrictions
remain: warm species need nonzero kpar; scans exclude kpar=0. No artificial
wavevector, temperature or frequency floors are used.

Corrected snapshots use the same real lab frequency `wr` for total and signed
harmonic susceptibilities. Centered susceptibility derivatives require stable
matrix and E-projected convergence; cold derivatives are analytic. The
accepted complex-frequency dispersion state and field/moment eigenfunctions
are preserved. Failed derivatives do not invalidate the dispersion root.

For E=deltaE/deltaEx and kappa=k*d_ref, define

```
HE = |E|^2
HB = |kappa cross E/(a*wr)|^2,  a=vA_ref/c
Hs = E† d[wr*chi_H,s(wr)]/dwr E
H  = HE + HB + sum_s Hs
W  = abs(deltaEx)^2/(16*pi) * H
Ps = Im(E† chi_s(wr) E)/H
```

Physical wavevectors, species parameters and reference scales remain fixed
under differentiation. No species Doppler frequency replaces lab wr in its
prefactor. Both kinetic and cold species contribute to Hs. Auxiliary thermal
scales are coordinates, not additional particles or energy terms.

`Ps` is a signed absorption/emission diagnostic, approximately a species
contribution to `-gamma/wr` for weak isolated modes. All seven historical
channels retain their definitions/order: n=0 yy, yz, zy, zz, full n=0,
transverse n=+1, transverse n=-1. The first four sum to full n=0; the last
three need not sum to total Ps. Summing all seven double-counts n=0 and also
omits higher harmonics and some +/-1 longitudinal/cross contributions.
Cold species' kinetic Ps/channels are unavailable (NaN), not zero. Fully cold
systems have energy but no kinetic-heating total.

The `.energy` companion reports individual energy components, signed
fractions, the noncancelling electric contribution `f_E_abs`, derivative
convergence, and numerical/applicability statuses. Negative H is allowed;
vanishing or numerically unresolved H has unavailable signed fractions.
These are formal real-frequency wave-energy expressions, not exact conserved
energies for arbitrary complex modes. `abs(gamma/wr)>0.05` and independently
failed weak-mode balance set applicability bits; normalized Ps/channels then
remain NaN, while resolved formal energy, raw Q, fields and moments remain
available. A weak-mode balance check alone cannot diagnose every degeneracy.

Legacy replay uses `heating_convention='legacy', z_backend='legacy',
bessel_mode='fixed128'`. It retains the historical arithmetic, normalizer,
shifted harmonic snapshots and collisional behavior. Explicit legacy heating
is unavailable in cold-containing scans; the direct evaluator can still
report fields/moments with its existing unavailable-heating message.
Corrected collisional heating is rejected rather than reinterpreting the
collision operator. Maps gain no heating, radial/plasma-parameter heating
remains legacy-only, and FPC behavior is unchanged.

Python `set_params(..., heating_convention='corrected')` writes and reads this
control. Read the companion with `load_plume_energy(filename)`; use the
existing `load_plume_sweep(..., nspec=N, eigen=..., heating=True)` for the
unchanged low_n=true columns. Supply switches explicitly for ambiguous files.

### Python wrapper roundtrips and unsupported requests

`plume_input.load_from_file(path)` preserves supported scalar controls in
`params`, `species_N`, `maps`, `scan_input_N`, `guess_N`, and `fpc`. This includes
`low_n`, both log flags, `nr`, `ni`, `tensor`, and all six Cartesian FPC bounds.
Loading replaces the previous object state. An explicit input `dataName` is
preserved; assign `obj.dataname` after loading to redirect output. Quoted paths,
inline scalar assignments, logical values and D exponents are supported.
Unknown controls, unsupported groups (including radial/k_range groups), and
array/expression syntax raise `ValueError` instead of being silently lost.
The wrapper is not a general Fortran namelist parser.

`set_maps(loggridw, omi, omf, gami, gamf, positive_roots,
loggridg=False, nr=128, ni=128)` keeps its historical first six arguments and
adds independent optional settings. `make_scan(..., tensor=False)` preserves
the tensor switch.

Exact-parallel **all-kinetic** runs (`kperp=0`) retain their historical
unsupported behavior and now warn once. No perpendicular floor is applied.
Mixed/cold analytic parallel limits remain available. Frequency maps do not
write eigenfunction or tensor diagnostics: requested flags now warn rather
than silently disappearing. Legacy map heating likewise warns and remains
absent. Corrected map heating retains its existing rejection. These warnings
preserve accepted legacy workflows; successful process exit alone does not
establish a valid unsupported calculation.
