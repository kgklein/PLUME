<a id="examples-a-twelve-lesson-curriculum"></a>

# Examples: a twelve-lesson curriculum

[Manual home](index.md) · [Quick Start](tutorial.md) · [Modes](modes.md)

All commands below assume the repository root, a completed `make -j2`, and
`mkdir -p data/example`. Links point to existing input sources; no new input
files are invented by this manual. Match each reader to the file's species count
and flags. Benchmarks below are current locally verified printed values or
explicitly labeled developer references, not universal dispersion relations.
A process finishing is not by itself validation of every mode/diagnostic.

<a id="lessons-with-existing-bounded-inputs"></a>
<a id="planned-lessons-do-not-present-missing-inputs-as-working-files"></a>

## Lessons with existing inputs

| Lesson | Input and physical purpose | Normalization / populations | Run and inspect |
|---|---|---|---|
| 1. Standard kinetic plasma and equivalent units | [normalization_rho_ref.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/normalization_rho_ref.in), [normalization_d_ref.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/normalization_d_ref.in) | Equivalent rho_ref/d_ref proton-electron inputs; alpha_R=1.5, electron drift .02 | Run both files; compare converted k columns, all three frequencies and E/B/moments; first w=.06024424-.0004186511i |
| 2. Low-frequency Alfvénic/kinetic-Alfvén response and channels | [heating_corrected_kinetic.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_kinetic.in) | rho_ref, two warm populations, anisotropic reference/electrons | Run file; inspect .mode1 and .energy; first w=.06098715-1.270997e-5i, H about 20342.68; identify branch by polarization, not filename alone |
| 3. Cold parallel ion-cyclotron/circular response | [cold_fields_parallel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_fields_parallel.in), [cold_proton_electron.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_proton_electron.in) | d_ref, stationary cold proton/electron, finite electron mass | Field scan: first w about .05820986 with roundoff imaginary part; compare cold circular n^2=R/L and transverse polarization. The single supplied-root file only prints roots |
| 4. Fully cold perpendicular electromagnetic mode | [cold_proton_electron_perpendicular.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/cold_proton_electron_perpendicular.in) | d_ref, fully cold, exact kpar=0, single wavevector | Run file; inspect printed root against O/X relations. Do not infer perpendicular scans or an Ex=1 O-mode eigenfunction are supported |
| 5. Warm ions / cold electrons | [heating_corrected_warm_reference.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_warm_reference.in) | d_ref, warm reference, cold electrons | .mode1 and .energy; first w=.06153842-5.564326e-6i, H about 20384.80; electron Ps/channels NaN while H_e is present |
| 6. Cold ions / warm electrons | [heating_corrected_cold_reference.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_cold_reference.in) | d_ref, betap=0, cold protons, warm electrons | First w=.05916735-1.588769e-6i, H about 20516.88; proton kinetic heat NaN, electron heat available; auxiliary scale is not physical temperature |
| 7. Several populations and parallel drift currents | [mixed_moments_multispecies.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/mixed_moments_multispecies.in) | d_ref, interleaved cold protons, warm electrons/alphas, cold alphas | First w=.05440212-1.366260e-6i; inspect every species' U/n; reconstruct current including density times equilibrium drift |
| 8. Temperature-anisotropy instability | [example_anisotropy_par.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_anisotropy_par.in) | `rho_ref`, warm protons/electrons; \(\beta_{\parallel,R}=3\), \(T_{\perp,p}/T_{\parallel,p}=3\), isotropic electrons, no drifts | Inspect `.mode1`, `.mode2` and `.energy`; locate intervals with \(\gamma>0\), identify the unstable branch by polarization, and check heating applicability before interpreting absorption |
| 9. Whistler/R branch | [example_anisotropy_par.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_anisotropy_par.in) | Same plasma and scan as lesson 8; two supplied frequency guesses | Follow the whistler/R branch using its frequency and field polarization as \(k_\parallel\) increases; check branch continuity and warnings rather than identifying the branch by mode slot alone |
| 10. Two-dimensional scale and beta scan | [example_simple_grid.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_simple_grid.in) | `rho_ref`, warm protons/electrons; logarithmic \(k_\perp\rho_R\) and \(\beta_{\parallel,R}\) axes at fixed \(k_\parallel\rho_R=10^{-3}\) | Inspect `simple_grid_kperp_betap.mode1`, `simple_grid.energy` and `simple_grid.warn`; reconstruct the 2D grid from coordinate columns and check energy status/applicability before interpreting heating |
| 11. Vacuum electric energy and fully cold energy | [heating_corrected_cold.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_corrected_cold.in) | d_ref, physical proton/electron mass ratio, stationary fully cold | First w about .05972252; H about 20607.63; H_E about 1.00721 and f_E about 4.88755e-5; kinetic heat remains unavailable |
| 12. Historical numerical/diagnostic replay | [heating_legacy.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_legacy.in) | rho_ref all kinetic; explicit legacy Z/fixed128/legacy heat | Inspect old mode Ps/channels, not a corrected .energy; record source revision and all flags |

For example:

```bash
./plume.e inputs/example/heating_corrected_warm_reference.in
./plume.e inputs/example/heating_corrected_cold_reference.in
./plume.e inputs/example/example_anisotropy_par.in
./plume.e inputs/example/example_simple_grid.in
```

Lessons 8 and 9 share a logarithmic scan from \(k_\parallel\rho_R=10^{-3}\)
to \(10^2\), at fixed \(k_\perp\rho_R=10^{-3}\). The input requests 400 output
steps and 32 substeps per output step for each of two branches, with eigenfunctions
and corrected heating enabled. This is a longer continuation example than the
bounded single-branch inputs above; inspect failed points and diagnostic statuses
throughout the scan.

Lesson 10 scans \(k_\perp\rho_R\) from \(10^{-3}\) to \(10^2\) and
\(\beta_{\parallel,R}\) from 1 to 10, with corrected heating and eigenfunctions.
Its 200 by 30 output intervals produce a \(201\times31\) grid, or 6,231 rows
per branch. The `.energy` companion records exact parameter states in comments.
Valid roots can have unavailable heating; inspect the companion statuses and
warning log throughout the grid.

Energy is only a formal weak-mode diagnostic; check energy_status and applicability
in the companion before interpreting normalized absorption. Lesson 11's ordinary
cold example has a small electric fraction. The **separate test-only** moderate-a,
artificial mass-ratio-4 case has f_E about .500000281; do not attribute that number
to the supplied physical-mass-ratio example. Creating its explicit public input
is later curated-example work.

The field-only warm/cold mixed files used by the [Quick Start](tutorial.md) are
simpler companions to lessons 5/6, not new physics examples. They are useful when
you want E/B/moments without introducing heating/energy interpretation.

<a id="secondary-numerical-and-historical-examples"></a>

## Secondary numerical and historical examples

- [adaptive_bessel.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/adaptive_bessel.in)
  and [faddeeva_accuracy.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/faddeeva_accuracy.in)
  are numerical stress/convergence exercises, not beginner benchmarks. Examine
  harmonic counts, cap warnings and tensor accuracy as described in [Numerics](numerics.md).
- [mixed_scan_warm_reference.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/mixed_scan_warm_reference.in)
  and [fully_cold_scan.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/fully_cold_scan.in)
  are existing serial continuation examples with no failed points in this replay.
- [mixed_scan_cold_reference.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/mixed_scan_cold_reference.in)
  finishes but has a failed written point in current replay. Study NaN/recovery
  deliberately; use mixed_fields_cold_reference.in for an all-success bounded
  first scan instead.
- [example_map_par.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_map_par.in)
  is an existing linear-map lesson; examine sampled determinants versus refined
  candidates and branch slots. It is larger than the Quick Start.
- [example_guess_double.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_guess_double.in)
  is a large historical 2D scan; the previous audit interrupted it at its short
  runtime limit. It is not certified as a small, all-success new-user example.
- Historical radial/parameter/FPC fixtures retain their specialized limitations.
  No new curated radial/FPC input or scientific extension is claimed here.

Further example work will automate benchmarks and document branch identity.
Existing inputs are not rewritten by this manual task. Current execution evidence and all linked-input existence
checks are recorded in the [manual report](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_5B1_MANUAL_REPORT.md).
