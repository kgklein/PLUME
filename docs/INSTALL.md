<a id="installation-and-verification"></a>

# Installation and verification

[Manual home](index.md) · [Next: Quick Start](tutorial.md)

<a id="requirements-by-activity"></a>

## Requirements by activity

| Activity | Requirements |
|---|---|
| Compile/run the solver | GNU Fortran, an ordinary C99 compiler, GNU make, standard Fortran/C runtime libraries and libm |
| Read/write inputs and analyze output w/ Python | Python 3, NumPy, Matplotlib |
| Field Particle Correlation analysis w/ Python | Python 3, NumPy, Matplotlib |
| Independent reference tests | SciPy, NumPy, Matplotlib; Python's unittest is built in |
| Optional high-precision tests | mpmath; tests fall back or skip the specific optional references when unavailable |
| Plotting | Matplotlib; historical `.plt` examples require gnuplot |
| Local documentation | FORD and its Python dependencies; Graphviz for API graphs |

The PLUME plasma dispersion function is bundled and compiled into the executable.
Users do not need SciPy, libcerf, a C++ runtime, CMake, Python, or a
package manager **to run the solver**.
Python was was used as an external, independent validation of some of the underlying methods.

GNU Fortran/C 11.4.0 and GNU make 4.3 have been exercised on Linux.
Recent validation used Python 3.10.12, NumPy 2.2.6, SciPy 1.15.3 and Matplotlib 3.10.8.
Other platforms/compilers are not certified by the presence of an alternative Makefile rule.
Preserve the project's promoted-double and strict-IEEE rules,
as changing precision or optimization can alter legacy regressions and validation.

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
built automatically. Do not link arbitrary kind/precision variants together.
The [numerical-method notes](numerics.md) explain the precision interface.

<a id="run-from-the-repository-root"></a>

## Run from the repository root

```bash
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

The command takes an input filename ending in `.in`.
Its indexed namelist extractor creates a scratch file beside that input,
so that directory must be writable.
It also has a historical 100-character argument buffer, so avoid very
long input paths. **This will be fixed in future revisions.**
The input's `dataName='example'` selects `data/example/`
relative to the **working directory**, not relative to the input file.

*PLUME does not create directories automatically.*
Files and numerical metadata may be opened even with `writeOut=.false.`
as this flag reduces `stdout` verbosity and is not a master file-output switch.
Use a unique `outputName` and an isolated directory for exploratory runs.
Files with matching names are replaced

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

The editable user-manual sources live in `docs/`, with the homepage in
`docs/index.md` and images in `docs/assets/`. The root `README.md` is the GitHub
project overview; `ford_project.md` contains the FORD configuration and API
landing-page body. Edit the Markdown sources rather than generated HTML.

Install FORD in a separate environment if desired (`python3 -m pip install ford`)
and provide Graphviz `dot` for graphs. From the repository root, build and open
the complete local manual with:

```bash
./scripts/build_docs.sh
xdg-open docs_out/page/index.html
```

The builder stages copies in a temporary directory, converts local Markdown
links, copies required assets, and validates the manual before replacing
`docs_out/`. The previous output is preserved in a uniquely named backup.
The same builder runs in documentation CI; building locally does not publish a
website. See [Developer documentation](development.md#local-documentation-build)
for validation details and remaining FORD limitations.
