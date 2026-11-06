<a id="installation-and-verification"></a>

# Installation and verification

[Manual home](./) · [Next: Quick Start](tutorial.md)

<a id="requirements-by-activity"></a>

## Requirements by activity

| Activity | Requirements |
|---|---|
| Compile/run the solver | GNU Fortran, an ordinary C99 compiler, GNU make, standard Fortran/C runtime libraries and libm |
| Read/write inputs and analyze output with Python | Python 3, NumPy; the current `linfpclib.linfpc` module also imports Matplotlib |
| Independent reference tests and some analysis | SciPy, NumPy, Matplotlib; Python's unittest is built in |
| Optional high-precision tests | mpmath; tests fall back or skip the specific optional references when unavailable |
| Plotting | Matplotlib; historical `.plt` examples require gnuplot; other plotting folders have their own assumptions |
| Local documentation | FORD and its Python dependencies; Graphviz for API graphs |

The accurate plasma dispersion function is bundled and compiled into the
executable. Users do not need SciPy, libcerf, a C++ runtime, CMake, Python, or a
package manager **to run the solver**. Python requirements apply to analysis,
wrapper and test workflows, not to susceptibility evaluation at runtime.

GNU Fortran/C 11.4.0 and GNU make 4.3 have been exercised on Linux. Recent
validation used Python 3.10.12, NumPy 2.2.6, SciPy 1.15.3 and Matplotlib 3.10.8.
Other platforms/compilers are not certified by the presence of an alternative
Makefile rule. Preserve the project's promoted-double and strict-IEEE rules:
changing precision or optimization can alter legacy regressions and validation.

<a id="get-the-source-and-compile"></a>

## Get the source and compile

```bash
git clone https://github.com/kgklein/PLUME.git
cd PLUME
make -j2
```

The resulting executable is `plume.e`. To rebuild from source:

```bash
make clean
make -j2
```

The default Makefile selects its GNU configuration. The existing flags promote
Fortran default real to binary64; several numerical modules are built with
strict floating-point semantics. The C99 Faddeeva sources and ABI bridge are
built automatically. Do not link arbitrary kind/precision variants together;
[numerical-method notes](numerics.md) explain the precision interface.

<a id="run-from-the-repository-root"></a>

## Run from the repository root

```bash
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

The command takes an input filename ending in `.in`. Its indexed namelist
extractor creates a scratch file beside that input, so that directory must be
writable. It also has a historical 100-character argument buffer: avoid very
long input paths. The input's `dataName='example'` selects `data/example/`
relative to the **working directory**, not relative to the input file.

PLUME does not create these directories automatically. Files and numerical
metadata may be opened even with `writeOut=.false.`. This flag reduces stdout
verbosity; it is not a master file-output switch. Use a unique `outputName` and
an isolated directory for exploratory runs: files with matching names can be
replaced, and old branch files can otherwise be mistaken for new output.

Follow the [Quick Start](tutorial.md) before attempting large maps or scans.

<a id="python-environment-and-tests"></a>

## Python environment and tests

In a virtual environment, install the analysis/test dependencies if needed:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install numpy scipy matplotlib
# Optional high-precision references:
python3 -m pip install mpmath
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

`MPLCONFIGDIR` is useful on systems whose home Matplotlib cache is not writable.
Tests compile additional isolated Fortran probes, including bounds-checked and
floating-point-trapping builds, and check pinned legacy revisions. Run in a Git
checkout with the historical commits available, not only an exported source
archive, for those comparisons. Optional arbitrary precision is a test reference,
not a production dependency. A current test count is recorded in the
[manual verification report](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_5B1_MANUAL_REPORT.md);
it is not a permanent scientific guarantee.

<a id="documentation-sources-and-local-builds"></a>

## Documentation sources and local builds

Root Markdown files are the editable user-manual sources. `ford_project.md`
contains the FORD project configuration and landing-page body. The existing
`.github/workflows/doc.yml` stages selected pages into `docs/` with `title:`
metadata, then runs:

```bash
ford ford_project.md
```

Install FORD in a separate environment if desired (`python3 -m pip install ford`)
and provide Graphviz `dot` for graphs. Running FORD directly in an unstaged
checkout does not stage the static manual pages for you. The current workflow
has a fixed staging list; several new manual pages are not yet included there.
Its link rewriting and assets also need follow-up. See
[Developer documentation](development.md#local-documentation-build) and the
manual verification report for exact locally tested builds and known failures.
Do not infer that a successful API build means the complete manual was included.

This manual update does not deploy a site or repair deployment configuration.
