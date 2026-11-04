# PLUME Input

This is a reference for the key input parameters used by PLUME.

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
`nscan=0`. Additional species may have exactly `tauS=0`; every positive
temperature, however small, remains kinetic. Cold `alphS` is ignored (zero is
allowed). Density, signed charge, mass and parallel drift retain their definitions.
No cold temperature is replaced by a numerical floor.
Parallel (`kperp=0`) and oblique propagation are supported. Mixed runs require
nonzero `kpar`: the warm susceptibility's exact perpendicular limit is not
implemented. No small nonzero wavevector is substituted automatically.

- `option=-1` evaluates the determinant and reports its status. For mixed runs
  it explicitly omits the usual supplementary eigenfunction/heating output.
- `option=1`, `use_map=.false.` refines supplied `&guess_N` values with the
  existing secant method. Set `writeOut=.true.` to print the solutions.
- `option=0` searches a frequency map. Singular samples remain NaN in map
  output and are excluded from root selection; failed refinements are skipped.
- Cold-containing parameter/radial scans, eigenfunctions, heating, and FPC
  are unavailable. Requests fail clearly; temperature scans that reach or leave
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
with `nscan=0`. The existing secant and map search are reused. Root guesses
select branches and must avoid exact poles; map resolution must resolve the
branch being sought. Eigenfunctions, heating (including harmonic-resolved
heating), FPC, collisions, and parameter/radial scans remain unavailable.
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
`collision_type=0`, and `nscan=0`. Additional species use exactly zero `tauS`
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

Cold-containing eigenfunctions, heating, harmonic heating, FPC, collisions,
parameter/radial scans and model-transition scans remain unavailable. Charge
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
Write or suppress output to screen.

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
Spacing automatically calculated between `omi` and `omf`.  

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{p}$ map search.
Spacing automatically calculated between `gami` and `gamf`  

**`omi`**  
Smallest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`omf`**  
Largest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_{p}$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_{p}$ value for complex map search.

**`nr`**  
Number of $\omega_{r}/\Omega_{p}$ points in frequency grid.
Default is 128 if not specified.

**`ni`**  
Number of $\gamma/\Omega_{p}$ points in frequency grid.
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
Calculates eigenfunctions if true.

**`tensor`**  
Outputs susceptibility tensor if true.     
