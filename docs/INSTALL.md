<a id="installation-and-verification"></a>

# Installation and verification

[Manual home](index.md) · [Next: Quick Start](tutorial.md)

<a id="requirements-by-activity"></a>

## System Requirements

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
Python packages were was used as an external, independent validation of some of the underlying methods.

*Caveats:* GNU Fortran/C 11.4.0 and GNU make 4.3 have been exercised on Linux.
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

Follow the [Quick Start](tutorial.md) before attempting anything major.

<a id="python-environment-and-tests"></a>

<a id="documentation-sources-and-local-builds"></a>

## Documentation sources and local builds

The editable user-manual sources live in `docs/`, with the homepage in
`docs/index.md` and images in `docs/assets/`. 
The root `README.md` is the GitHub project overview, 
and `ford_project.md` contains the FORD configuration and API landing-page body. 
Edit the Markdown sources, not the generated HTML.

Install FORD in a separate environment if desired (`python3 -m pip install ford`)
and provide Graphviz `dot` for graphs. 
From the repository root, build and open the complete local manual with:

```bash
./scripts/build_docs.sh
xdg-open docs_out/page/index.html
```

The builder stages copies in a temporary directory, converts local Markdown
links, copies required assets, and validates the manual before replacing
`docs_out/`. 
The previous output is preserved in a uniquely named backup.
The same builder runs in documentation CI; building locally does not publish a
website. See [Developer documentation](development.md#local-documentation-build)
for validation details and remaining FORD limitations.
