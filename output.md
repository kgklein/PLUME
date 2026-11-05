# PLUME Output

## Cold-reference mixed output (Stage 3B.2)

Supplied roots and frequency maps retain `omega/Omega_ref` and native `d_ref`
input units. The direct evaluator prints raw `det(Lambda)`, as for fully cold
runs. The normalization sidecar also records `auxiliary_first_warm_species`
and `b_star`, selected from its physical temperature; the exact-parallel limit
uses physical species thermal speeds directly. This metadata is not a change
to any physical input or a new user parameter.

`susc` stores the ordinary response for each species and `lam` stores the
complete wave matrix. Status 5 denotes unrepresentable/nonfinite arithmetic,
including auxiliary conversion or kinetic intermediates. Failure invalidates
all response arrays; isolated failed map samples stay NaN and are excluded from
candidate selection. No heating or eigenfunction values are fabricated.
See the [cold-reference report](tests/COLD_REFERENCE_MIXED_REPORT.md).

## Fully cold output (Stage 3A)

Fully cold `d_ref` runs use the same root/map layouts, frequency units,
status codes and failure recovery described below. The direct evaluator
(`option=-1`) prints **raw `det(Lambda)`**, explicitly labeled, because the
legacy thermal output scale is undefined at zero reference temperature.
The normalization sidecar records this distinction and the absence of thermal
reference scales. The Fortran `disp` API always returns the raw determinant.
Cold-containing diagnostic limitations apply to fully cold runs too.
See the [Stage 3A report](tests/FULLY_COLD_DISPERSION_REPORT.md).

## Mixed kinetic/cold output (Stage 2)

Mixed `d_ref` runs retain the frequency normalization `omega/Omega_R` and
ordinary species susceptibility `chi_s`. Supplied-guess runs (`option=1`,
`nscan=0`) print roots when `writeOut=.true.`. Frequency maps retain their numeric
layout and normalization sidecar. Singular/nonfinite map samples have NaN
determinant and log-magnitude fields, are excluded from root selection, and
produce a summary warning. Failed candidate refinements are omitted.

Direct evaluation (`option=-1`) retains the legacy printed determinant scaling
`det(Lambda)*solver_vtp**6` and adds a mixed status line: 0 valid, 1 invalid input,
2 zero frequency, 3 zero Doppler frequency, 4 cyclotron pole, 5 nonfinite numerical
result. The Fortran `disp` API returns the unscaled determinant; inspect
`vars::dispersion_status` before using its tensor arrays after a mixed call.
Failure invalidates `susc`, `lam`, and resonance arrays rather than leaving
stale values. A subsequent valid evaluation recomputes the physical response.

Mixed eigenfunctions, heating and FPC are unavailable; the direct evaluator
prints an explicit omission notice, and dedicated diagnostic requests fail.
Cold resonance-resolved susceptibility entries (`susc_low`) are undefined/NaN,
not fictitious heating contributions. See the [Stage 2 report](tests/MIXED_DISPERSION_REPORT.md).

PLUME writes output solutions to the `/data` directory, to the subdirectory specified by `dataName`.
All output file names will have the `outputName` string included in the name to distinguish between distinct calculations. 

Numeric output columns contain values in the convention selected by
the input `normalization` parameter.
The first four columns of roots, modes, and tensor files are
`kperp`, `kpar`, `betap`, and `vtp` in **input units**:

- $k\rho_R$ and $v_{t\parallel,R}/c$ for `rho_ref`,
- or $kd_R$ and $v_{A,R}/c$ for `d_ref`.

Species `tauS` columns likewise contain the selected temperature definition.

The complex frequency always uses $\Omega_R$.
A `d_ref` run also writes `normalization_<outputName>.txt` beside its numeric outputs;
its absence in older outputs indicates the legacy `rho_ref` convention.
Reusing an output name with `rho_ref` removes a stale `d_ref` metadata file.
Existing numeric file column order and formatting are unchanged.

For the gyrotropic and Cartesian FPC files, velocities and fields continue to
use the kinetic solver's thermal velocity convention;
their headers describe those derived quantities.
The normalization metadata identifies the input convention used to generate the run.

## dispersion_*outputName*.map

Value of the dispersion tensor $\mathcal{\Lambda}(\omega_{\textrm{r}},\gamma)$ on a defined complex frequency grid.  
Solutions to the dispersion relation satisfy $|\mathcal{\Lambda}|  =0$.
This file is generated from the *map_search* subroutine in *disprels.f90*,
and invoked when `use_map` =.true. 

The data is ordered in columns as  
1. $\omega_r$  
2. $\gamma$   
3. $\log_{10} |\mathcal{\Lambda}|$  
4. Re $[|\mathcal{\Lambda}|]$  
5. Im $[|\mathcal{\Lambda}|]$  

The *&maps* namelist in *filename*.in determines the structure of *filename*.map.  
The range of $\omega_{\textrm{r}}/\Omega_p$ is from `omi` to `omf` with `nr` steps.
Logorithmic or linear spacing is selected with `loggridw`.
The range of $\gamma_{\textrm{r}}/\Omega_p$ is from `gami` to `gamf` with `ni` steps. Logorithmic or linear spacing is selected with `loggridg`.
`nr` and `ni` are set to 128 as a default which can be replaced by user defined parameters;
see [input.md].

## dispersion_*outputName*.roots

Identified solutions to the dispersion relation $|\mathcal{D}|  =0$, calculated using *refine_guess* in *disprels.f90*.

The data is ordered in columns as:
1. $k_\perp \rho_{ref}$
2. $k_\parallel \rho_{ref}$
3. $\beta_{ref,\parallel}$
4. $v_{t,ref,\parallel}/c$
5. $\omega_{r}/\Omega_{ref}$
6. $\gamma/\Omega_{ref}$

This will be followed by 6`nspec` columns containing the parameter lists $\mathcal{P}_j$ for each species or component.
- $T_{ref,\parallel}/T_{j,\parallel}$.
- $m_{ref}/m_{j}$.
- $T_{j,\perp}/T_{j,\parallel}$.
- $q_{ref}/q_{j}$.
- $n_{j}/n_{ref}$.
- $v_{j,drift}/v_{A,ref}$.

Only the first `nroot_max` solutions will be identified and written to file.

## *outputName*_*param*_*range*.modeN

The dispersion relation calculated from `om_scan` will have different formats based upon the number of species as well as the suplementary calculations invoked.

The *param* name will be assigned based on the kind of parameter scan performed, specifically:

*`scan_style`=-1:*

- `k` for scan from $\textbf{k}_0 \rho_{ref}$ to $\textbf{k}_1 \rho_{ref}$ (`scan_type`=0).
- `theta` for scan from $\theta_0$ to $\theta_1$ with fixed $|k|\rho_{ref}$ (`scan_type`=1).
- `k_fixed_theta` for scan from $|k|_1\rho_{ref}$ to $|k|_1\rho_{ref}$ with a fixed value of $\theta$ (`scan_type`=2).

*`scan_style`=0:*

- `kperp` for scan of $k_\perp \rho_{ref}$ (`scan_type`=0).
- `kpar` for scan of $k_\parallel \rho_{ref}$ (`scan_type`=1).
- `beta` for scan of $\beta_{ref,\parallel}$ (`scan_type`=2).
- `vtp` for scan of $v_{t,ref,\parallel}/c$ (`scan_type`=3).

*`scan_style`=1-`nspec`:*

- `tauS` for scan of $T_{ref,\parallel}/T_{j,\parallel}$ (`scan_type`=0).
- `muS` for scan of $m_{ref}/m_{j}$ (`scan_type`=1).
- `alphS` for scan of $T_{j,\perp}/T_{j,\parallel}$ (`scan_type`=2).
- `qs` for scan of $q_{ref}/q_{j}$ (`scan_type`=3).
- `Ds` for scan of $n_{j}/n_{ref}$ (`scan_type`=4).
- `Vs` for scan of $v_{j,drift}/v_{A,ref}$ (`scan_type`=5).

The initial and final values of the parameter scanned over will be reflected in the *range* portion of the file name.

One file for each mode followed will be produced.

The first six columns are the same regardless of supplemental calculations,
normalized to the appropriate convention.
They are ordered as
1. $k_\perp L_{ref}$
2. $k_\parallel L_{ref}$
3. $\beta_{ref,\parallel}$
4. $v_{t,ref,\parallel}/c$ or $v_{A,ref}/c$
5. $\omega_{\textrm{r}}/\Omega_{ref}$   
6. $\gamma/\Omega_{ref}$

If `eigen` is set to true, the next set of columns will be the eigenfluctuations, with fields normalized to $E_x$, velocity normalized to $c E_x/B_0$ and density normalized to $n_{0s} E_x/B_0$, 

7. Re $[B_x]$   
8. Im $[B_x]$   
9. Re $[B_y]$   
10. Im $[B_y]$   
11. Re $[B_z]$   
12. Im $[B_z]$   
13. Re $[E_x]$   
14. Im $[E_x]$   
15. Re $[E_y]$   
16. Im $[E_y]$   
17. Re $[E_z]$   
18. Im $[E_z]$   
19. [+6(j-1)] Re $[\delta U_{x,j}]$   
20. [+6(j-1)] Im $[\delta U_{x,j}]$   
21. [+6(j-1)] Re $[\delta U_{y,j}]$   
22. [+6(j-1)] Im $[\delta U_{y,j}]$   
23. [+6(j-1)] Re $[\delta U_{z,j}]$   
24. [+6(j-1)] Im $[\delta U_{z,j}]$   
19. [+6(`nspec`)+2(j-1)] Re $[\delta n_{j}]$   
20. [+6(`nspec`)+2(j-1)] Im $[\delta n_{j}]$ 

where `j` ranges from 1 to `nspec`.
The normalization follows Eqns. 33-37 in [Klein, K. G., Howes, G. G.,
and Brown, C. R., 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)

If `heat` is set to true, the next set of columns will be the power absorption or emission from each component. If `low_n` is set to true, additional terms associated with Landau, Transit time, and Cyclotron heating will be output. If `eigen` is false, this data will start in the 7th column. If `eigen` is true, this data will start in the 18+8 `nspec`+1st column.

If `low_n` is true, the ordering is:
  1) All total powers P_j for j=1..nspec, then
  2) For each species j, the six sub-terms (yy, yz, zy, zz, n=0, n=±1) consecutively.

Concretely (1-based column numbers; {!eigen}=1 if eigen is false, else 0):

- P_j: 18 + (8 - 6{!eigen})*nspec + j

Let afterP = 18 + (8 - 6{!eigen})*nspec + nspec.
Then for species j in {1,..,nspec}:

- P_j^yy:   afterP + 6*(j-1) + 1   (Transit Time Damping term 1)
- P_j^yz:   afterP + 6*(j-1) + 2   (Transit Time Damping term 2)
- P_j^zy:   afterP + 6*(j-1) + 3   (Landau Damping term 1)
- P_j^zz:   afterP + 6*(j-1) + 4   (Landau Damping term 2)
- P_j^{n=0}: afterP + 6*(j-1) + 5  (sum of Landau and Transit Time Damping)
- P_j^{n=+1}: afterP + 6*(j-1) + 6 (n=±1 Cyclotron Damping)
- P_j^{n=-1}: afterP + 7*(j-1) + 7 (n=-1 Cyclotron Damping)

Here, {!eigen} (negation of the eigen boolean) equals 1 if eigen is false and 0 if eigen is true.


This will be followed by 6`nspec` columns containing the parameter lists $\mathcal{P}_j$ for each species or component. 
- 19+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $T_{ref,\parallel}/T_{j,\parallel}$.
- 20+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $m_{ref}/m_{j}$.
- 21+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $T_{j,\perp}/T_{j,\parallel}$.
- 22+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $q_{ref}/q_{j}$.
- 23+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $n_{j}/n_{ref}$.
- 24+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $v_{j,drift}/v_{A,ref}$.
Here, noutperspec = 0 is the number of additional outputs created by setting heating or eigen to true. If eigen and heating are false, then noutperspec = 8, if eigen is false and heating is true, then noutperspec = 8, and if eigen is true and if eigen and heating are true, then noutperspec=16. Note that {!eigen} (negation of eigen boolean) is equal to 1 if eigen is false and 0 if eigen is true.

For ease of reference, consider the below tables for column numbers for two and three component calculations.
   
![Table for output column numbering for a two-component plasma.](./data/PLUME-Output-Table.svg)

![Table for output column numbering for a three-component plasma.](./data/PLUME-Output-Table-2.svg)

![Table for output column numbering for a four-component plasma.](./data/PLUME-Output-Table-3.svg)

This same data structure is preserved for the output from `om_double_scan`.
The file naming convention will include a value of *param* from both of the two parameters scanned, and the code will not output information about the range of parameters in the file name.


## Bessel harmonic metadata

Normal executable completion writes
`data/<dataName>/<outputName>.bessel`, without changing numeric mode, root,
tensor, heating or FPC file layouts. It records harmonic mode, tolerance,
adaptive maximum, fixed count, and whether unresolved truncation was encountered.
With `bessel_diagnostics=.true.`, or whenever an unresolved adaptive sum occurred,
it includes:

- Species index and physical Bessel argument.
- Active maximum absolute harmonic order and allocated capacity (which can be larger).
- Convergence status and upper bounds for ordinary and quadratic-weighted tails.
- Cache hits, misses and coefficient-array evaluations, counted per species slot.
- History of capped adaptive requests, including those preceding a later resolved point.

A miss can require selection only: a changed tolerance may reuse previously
computed coefficients. An array evaluation calculates all required Bessel
coefficients together, not one coefficient. Density, drift, parallel wavevector,
frequency and auxiliary scale do not enter the physical cache key. Exact argument
matching is used; no approximate matching tolerance is applied. A species-count
change preserves surviving slots, while a different argument in a reordered slot
requires recalculation. No on-disk cache is read.

Fixed 128 preserves its historical evaluators and caches; detailed hit/miss
instrumentation is unavailable for that compatibility path and is labeled as
such. Its final-point tail bounds still determine the `unresolved` flag.
Analytical parallel and fully cold paths allocate no adaptive Bessel arrays.
`unresolved=F` describes truncation only, not a bound on total numerical error.
For standalone library/probe callers, `finish_bessel_diagnostics()` writes the
sidecar and `bessel_cache::write_bessel_stats(unit)` provides the table directly.


## Collisional longitudinal diagnostics

The collisional `zz` susceptibility and its `susc_low` entries for n=0,+1,-1
now use the corrected longitudinal Krook response. Collisional roots, fields
and heating derived from that response can change. The correction is independent of the selected Z backend.
The signed `zz` diagnostics contain each harmonic's full nonresonant terms;
the bulk sum separately regroups the drift contribution. For a truncated sum,
reconstruction must include the drift remainder proportional to
`1-Gamma_0-2*sum(Gamma_n)`, rather than adding another full drift term.
See the [derivation and validation report](tests/COLLISIONAL_ZZ_CORRECTION.md).


## Plasma dispersion backend metadata and failures

With `writeOut=.true.`, input validation prints `Z backend: faddeeva` or
`Z backend: legacy` once. Bessel sidecars append `# z_backend=<mode>` without
changing the existing numerical table rows. No new columns are added to numeric
root, tensor, eigenfunction or heating output files.

The accurate backend returns invalid determinants/tensors on nonfinite or
out-of-range evaluation, with the existing dispersion failure status. Failed
root trials are rejected; isolated invalid map points are excluded from minima
searches and retained as invalid values in map output. A subsequent valid trial
can succeed normally. Direct single-frequency failures skip supplementary
diagnostics. Compatibility mode preserves historical Z failure behavior.
See [backend selection](input.md#plasma-dispersion-function-backend) and the
[validation report](tests/FADDEEVA_BACKEND_REPORT.md).
