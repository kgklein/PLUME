<a id="plume-plasma-in-a-linear-uniform-magnetized-environment"></a>

# PLUME: Plasma in a Linear Uniform Magnetized Environment

![PLUME logo](./docs/assets/PLUME_logo.png)

PLUME is a serial Fortran solver for linear electromagnetic dispersion relations in homogeneous, nonrelativistic magnetized plasmas.
Populations may be combinations of drifting kinetic bi-Maxwellians or cold fluids.
Parallel drifts and kinetic temperature anisotropy are supported.

Given a wavevector and plasma state, PLUME calculates complex mode frequencies,
continues solutions through wavevector and/or parameter scans, and outputs electromagnetic
polarizations, species density/velocity perturbations, and signed absorption/emission for kinetic populations.
Details on
[models](docs/models.md) and [limitations](docs/limitations.md)
are available.

<a id="start-here"></a>

## Start here

From the repository root:

```bash
make -j2
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

This ten-output-point scan uses warm protons and cold electrons.
It writes a mode file, 'data/example/mixed_fields_warm_reference_kpar_60_90.mode1'
containing the complex frequency, electromagnetic polarization and species moments.
The [Quick Start](docs/tutorial.md) explains the input, expected result and output columns.
*The input directory must be writable and output directories must exist,
as `writeOut` controls verbosity, not file creation.*

<a id="user-manual"></a>

## Documentation

The [user manual](docs/index.md) provides the complete contents and scientific references.
Begin with [installation](docs/INSTALL.md) and the [Quick Start](docs/tutorial.md).
More details can be found in the [input reference](docs/input.md),
[normalization](docs/normalization.md), [calculation modes](docs/modes.md) and
[output guide](docs/output.md) sections for specific calculations.
The [examples](docs/examples.md) provide a longer introduction.
The [reproducibility guide](docs/reproducibility.md) documents historical numerical choices from PLUME 1.y.z. 

<a id="scientific-conventions"></a>

## Scientific conventions

PLUME uses Gaussian electromagnetic units and perturbations proportional to
\(\exp(i\mathbf k\cdot\mathbf r-i\omega t)\).
Solutoins with positive \(\operatorname{Im}\omega\) denotes growth.
Input wavevectors and temperatures follow either `rho_ref` or `d_ref` normalization.
The `d_ref` normalization is necessary to use cold species contributions; see
[normalization](docs/normalization.md) for details.

<a id="authors-scientific-references-and-citation"></a>

## Citation and license

For publications, cite the [PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2),
acknowledge the [repository](https://github.com/kgklein/PLUME), and record the
source revision and numerical settings used. The manual contains the
[scientific references and BibTeX entry](docs/index.md#authors-scientific-references-and-citation).

PLUME is distributed under the [BSD 2-Clause license](LICENSE),
with the bundled Faddeeva implementation retaining its separate [MIT license](third_party/faddeeva/README.md).
