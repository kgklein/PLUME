<a id="plume-plasma-in-a-linear-uniform-magnetized-environment"></a>

# PLUME: Plasma in a Linear Uniform Magnetized Environment

![PLUME logo](./docs/assets/PLUME_logo.png)

PLUME is a serial Fortran solver for linear electromagnetic dispersion relations
in homogeneous, nonrelativistic magnetized plasmas. Populations may be drifting
kinetic bi-Maxwellians or collisionless cold fluids, including mixtures of the
two. Parallel drifts and kinetic temperature anisotropy are supported; the
physical reference population is always `species_1`.

Given a wavevector and plasma state, PLUME calculates complex mode frequencies,
continues solutions through wavevector scans, and outputs electromagnetic
polarization and species density/velocity perturbations. Additional diagnostics
include signed absorption/emission for kinetic populations, with seven selected
harmonic/energy channels, and formal wave energy with vacuum field and species
contributions. Not every diagnostic is available for every plasma model; see
[species and models](docs/models.md) and [limitations](docs/limitations.md).

<a id="start-here"></a>

## Start here

From the repository root:

```bash
make -j2
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

This three-output-point scan uses warm protons and cold electrons. It writes a
mode file containing the complex frequency, electromagnetic polarization and
species moments. The [Quick Start](docs/tutorial.md) explains the input,
expected result and output columns. The input directory must be writable and
output directories must exist; `writeOut` controls verbosity, not file creation.

<a id="user-manual"></a>

## Documentation

The [user manual](docs/index.md) provides the complete contents and scientific
references. Begin with [installation](docs/INSTALL.md) and the
[Quick Start](docs/tutorial.md); consult the [input reference](docs/input.md),
[normalization](docs/normalization.md), [calculation modes](docs/modes.md) and
[output guide](docs/output.md) for individual calculations. The
[examples](docs/examples.md) provide a longer introduction, while the
[reproducibility guide](docs/reproducibility.md) documents historical numerical
choices. The [published manual](https://kgklein.github.io/PLUME/) may lag the
current checkout.

<a id="scientific-conventions"></a>

## Scientific conventions

PLUME uses Gaussian electromagnetic units and perturbations proportional to
\(\exp(i\mathbf k\cdot\mathbf r-i\omega t)\); positive
\(\operatorname{Im}\omega\) denotes growth. Input wavevectors and temperatures
follow either `rho_ref` or `d_ref` normalization. These conventions are not
interchangeable by changing only the normalization flag; see
[normalization](docs/normalization.md). PLUME does not support relativistic
kinetics, arbitrary distribution functions, or general collision operators;
exact-axis and diagnostic restrictions are described in the
[limitations](docs/limitations.md).

<a id="authors-scientific-references-and-citation"></a>

## Citation and license

For publications, cite the [PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2),
acknowledge the [repository](https://github.com/kgklein/PLUME), and record the
source revision and numerical settings used. The
[archived v1.0.1 release](https://doi.org/10.5281/zenodo.15215514) predates
some current capabilities. The manual contains the
[scientific references and BibTeX entry](docs/index.md#authors-scientific-references-and-citation).

PLUME is distributed under the [BSD 2-Clause license](LICENSE); the bundled
Faddeeva implementation has a separate [MIT license](third_party/faddeeva/README.md).
