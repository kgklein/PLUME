<a id="quick-start-a-small-mixed-plasma-wave"></a>

# Quick Start: a small mixed-plasma wave

[Manual home](index.md) · [Installation](INSTALL.md) · [Input reference](input.md)

For a first calculation, let's use the
[mixed_fields_warm_reference.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/mixed_fields_warm_reference.in).
This simulates a warm, isotropic proton VDF with cold electrons.
The solution is at an oblique wavevector, 
tracking one supplied complex-frequency guess for the solution, 
and follows along ten steps in wavevector space. 

<a id="1-build-and-run"></a>

## 1. Build and run

From the repository root:

```bash
make -j2
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

The input directory must be writable for the indexed-namelist scratch file.
The supplied example has  `writeOut=.true.`, so it also prints root information.
The primary file is:

```
data/example/mixed_fields_warm_reference_kpar_60_90.mode1
```

The integer filename tags encode scan endpoints.
Read the numeric columns in the file for details. 
Produced metadata include
`normalization_mixed_fields_warm_reference.txt` and
`mixed_fields_warm_reference.bessel`. 
The exact normalization metadata details are documented in the [output reference](output.md#metadata).

<a id="2-read-the-input-in-physical-terms"></a>

## 2. Read the input in physical terms

This compact rendering has the same essential state and scan as the existing
example; its optional Bessel diagnostics are omitted here for readability:

```fortran
&params
normalization='d_ref', betap=.4, vtp=.01
kperp=.05, kpar=.06
nspec=2, nscan=1, option=1, nroot_max=1
use_map=.false., writeOut=.true.
dataName='example', outputName='mixed_fields_warm_reference'
/
&species_1
 tauS=.4, muS=1, alphS=1.3, Qs=1, Ds=1, vvS=0
/
&species_2
 tauS=0, muS=1836, alphS=0, Qs=-1, Ds=1, vvS=0
/
&guess_1
 g_om=.06, g_gam=-.0001
/
&scan_input_1
 scan_style=0, scan_type=1, swi=.06, swf=.09
 ns=10, nres=2, swlog=.false.
 eigen=.true., heating=.false., tensor=.false.
/
```

Keep each indexed `&species_N`, `&guess_N` and `&scan_input_N` header **alone on
its own line**, in lowercase. The Fortran extractor is more restrictive than
ordinary case-insensitive namelist field syntax.

Here `vtp=.01` means \(a=v_{A,R}/c=.01\), and `kpar`/`kperp` mean
\(k_\parallel d_R\), \(k_\perp d_R\). The warm proton `tauS=.4` is its
Alfvénic parallel temperature and equals reference `betap`. Electron `tauS=0`
selects a true cold-fluid response; its anisotropy placeholder is ignored.
Signed `Qs=-1` gives negative electron gyrofrequency. 

`option=1` refines the guess and then continues a branch through serial scans.
`ns=10` means **ten output intervals, eleven written points**. 
`nres=2` inserts hidden continuation substeps. 
The written parallel coordinates are linear steps between .06 and .09; the physical plasma state and perpendicular scale remain fixed.
For how to vary other parameters, see [input reference](input.md)
The default plasma dispersion function Z backend is Faddeeva, 
with adaptive Bessel selection, 
and a corrected heating convention.

<a id="3-inspect-the-first-result"></a>

## 3. Inspect the first result

For the example input, the first row should give approximately

\[
 w=\frac{\omega}{\Omega_R}=0.06153842-5.564326\times10^{-6}i.
\]

The negative imaginary part signifies damping with time, not spatial attenuation.
The run should produce eleven rows with finite frequencies, no `scan root failed` warning, and finite E/B/moments.
The code exiting alone does not certify a valid mode or branch.

The file has 46 columns: six coordinates/root entries, 12 field entries,
16 moment entries and 12 species parameters. Columns 1–6 (one-based) are native
kperp, kpar, betap, vtp, real frequency, imaginary frequency. Columns 7–12 are
B real/imaginary pairs; 13–18 are E pairs. The remaining blocks are described
without hard-coded species assumptions in the [output schema](output.md#mode-files).

<a id="4-read-with-python"></a>

## 4. Read with Python

```python
import numpy as np
from linfpclib.linfpc import load_plume_sweep

path = 'data/example/mixed_fields_warm_reference_kpar_60_90.mode1'
d = load_plume_sweep(path, nspec=2, eigen=True, heating=False)
w = d['w'] + 1j*d['g']
print(w)
assert len(w) == 11 and np.isfinite(w).all()
print('Ex:', d['exr'][0] + 1j*d['exi'][0])  # 1+0j
```

The Ex=1 convention fixes complex amplitude and phase. E is divided by deltaEx,
B by the same deltaEx in Gaussian units. 
See [Eigenfunctions and moments](eigenfunctions.md) before
comparing these amplitudes with dimensional observations.

<a id="5-make-the-next-calculation-deliberately"></a>

## 5. Next Steps

- For all-kinetic absorption and wave energy, run
  `./plume.e inputs/example/heating_corrected_kinetic.in` and inspect its
  `.energy` companion, as explained in [Heating and wave energy](heating.md).
- For a cold physical reference, run
  `./plume.e inputs/example/mixed_fields_cold_reference.in`.
- For equivalent physical inputs in both normalizations, use the
  [matched pair](normalization.md#worked-equivalent-inputs).
- For maps or 2D scans, read [Calculation modes](modes.md) and use the
  [examples index](examples.md).

If no file appears, check the working directory, `dataName`, flags and input file.
For NaNs or missing eigenfunctions, inspect warnings and diagnostic statuses;
[troubleshooting](limitations.md) helps to distinguish a "bad root" from an unavailable
supplementary quantity.
