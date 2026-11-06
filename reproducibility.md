<a id="historical-reproduction-and-provenance"></a>

# Historical reproduction and provenance

[Manual home](./) · [Numerical methods](numerics.md) · [Heating](heating.md)

There are four different meanings of “legacy.” Record which one you need.
None of the input switches selects the others automatically.

| Choice | Setting | What is preserved |
|---|---|---|
| Z evaluator | `z_backend='legacy'` | Historical zetout numerical arithmetic/stopping/failure limitations |
| Harmonic sum | `bessel_mode='fixed128'`, `bessel_nfixed=128` | Historical maximum absolute order 128, i.e. -128..128, and original fixed-context arithmetic |
| Heating convention | `heating_convention='legacy'` | Historical all-kinetic normalizer, derivative, harmonic/scratch and output behavior |
| Physical/source history | A specific Git commit, compiler and flags | The equations, dispatch, bug fixes, formats and build rules in that revision |

Current defaults are Faddeeva/adaptive/corrected. Old files omitting the
selectors use **current defaults**, even though omitted normalization still
defaults to rho_ref. Selecting fixed128 alone does not reproduce historical Z
or heating. A legacy file should state all needed controls explicitly:

```fortran
&params
normalization='rho_ref'
z_backend='legacy'
bessel_mode='fixed128'
bessel_nfixed=128
heating_convention='legacy'
! Supply the remaining physical/workflow inputs as in the working file below.
/
```

That snippet is a selector fragment, not a complete runnable plasma. The existing
[heating_legacy.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/heating_legacy.in)
is complete and runnable:

```bash
make -j2
mkdir -p data/example
./plume.e inputs/example/heating_legacy.in
```

It produces the old all-kinetic mode heating blocks; do not expect corrected
vacuum-electric energy or a newly written .energy companion. Legacy cold/mixed
heating is unavailable. The old Ew slot is not a trustworthy wave energy.
Legacy Z intentionally retains approximation/failure defects, rather than being
an accurate fallback for an invalid modern evaluation.

<a id="what-numerical-selectors-do-not-undo"></a>

## What numerical selectors do not undo

These switches do not revert the separately corrected collisional zz physics,
new validation/warnings, new model dispatch, Python reader corrections, native
input definitions, or historical logarithmic-map coordinates. They do not
restore arbitrary old file formats or turn off all later behavior. For exact
old physics or byte-level reproduction of a different release, use its source
revision and compatible toolchain in a separate checkout. Do not overwrite a
working research tree to emulate an older calculation.

The existing exact tests pin pre-normalization commit
`c6c7616e1fc58eb6118e1f9423a3f0df6d52ef20` for two collisionless numerical
fixtures, and `8331925c7953dcc68f29f14c1d7b78dd9198bd5a` for kinetic scan/field
compatibility. Those are specific verified comparisons with explicit settings,
not a promise of universal byte identity for every input/compiler/platform.

<a id="keep-a-reproducibility-record"></a>

## Keep a reproducibility record

Save the full input, Git commit, compiler/version/flags, numerical selectors,
normalization, original species order, stdout warnings and output companions.
For scans record bounds, ns/nres, initial guess/map and branch identity evidence.
Retain unresolved harmonic/energy flags; do not strip failed rows to improve
plots. The normalization pair demonstrates equivalent physics after conversion,
not byte-identical native coordinate columns.

The versioned energy header and Bessel metadata help, but neither replaces the
input and source record. Most mode/root files have no self-describing header;
absence of normalization metadata does not identify arbitrary external files.
A source archive without Git history cannot run the pinned-revision tests as
if it were a full checkout.
