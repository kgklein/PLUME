<a id="plume-plasma-in-a-linear-uniform-magnetized-environment"></a>

# PLUME: Plasma in a Linear Uniform Magnetized Environment

![PLUME logo](./docs/assets/PLUME_logo.png)

PLUME calculates **linear electromagnetic dispersion relations in homogeneous,
nonrelativistic magnetized plasmas**. Each population is a drifting kinetic
bi-Maxwellian or a collisionless cold fluid. You can combine multiple populations,
including populations of the same particle species, with parallel drifts and
kinetic temperature anisotropy. The physical reference population is always
`species_1`.

Given a wavevector and plasma state, PLUME searches for complex frequencies,
continues modes through wavevector scans, and calculates electromagnetic
polarization and species density/velocity perturbations. Collisionless kinetic
populations also have signed absorption/emission diagnostics with seven selected
harmonic/energy channels. Corrected wave-energy diagnostics include vacuum
electric, magnetic and individual species contributions, including cold species.

PLUME is a **serial** Fortran executable. It does not provide relativistic
kinetics, arbitrary velocity distributions, general collision operators, or
universal exact-parallel/perpendicular support. The available Krook operator and
FPC pathways have additional restrictions. Start with the capability table in
[Species and models](docs/models.md) and the [limitations](docs/limitations.md).

<a id="start-here"></a>

## Start here

```bash
make -j2
mkdir -p data/example
./plume.e inputs/example/mixed_fields_warm_reference.in
```

This small three-output-point scan uses warm protons and cold electrons and
writes a mode file containing frequency, E/B and species moments. The
[Quick Start](docs/tutorial.md) explains its input, expected first root, columns and
success checks. Run from the repository root; the input directory must be
writable and output directories must exist. `writeOut` controls verbosity,
not whether files are produced.

<a id="user-manual"></a>

## User manual

| Your question | Authoritative page |
|---|---|
| How do I build and test PLUME? | [Installation](docs/INSTALL.md) |
| How do I run and read a first calculation? | [Quick Start](docs/tutorial.md) |
| What does every input control mean? | [Input reference](docs/input.md) |
| What are the physical units and reference scales? | [Normalization](docs/normalization.md) |
| Which kinetic, cold and mixed states are supported? | [Species and models](docs/models.md) |
| How do roots, maps and scans work? | [Calculation modes and scans](docs/modes.md) |
| How are fields and species moments normalized? | [Eigenfunctions and moments](docs/eigenfunctions.md) |
| What do absorption and formal wave energy mean? | [Heating and wave energy](docs/heating.md) |
| Which file/reader contains my diagnostic? | [Output files and Python](docs/output.md) |
| What should I study next? | [Twelve-lesson examples index](docs/examples.md) |
| What numerical methods and backends are used? | [Numerical methods](docs/numerics.md) |
| Why is a result missing or NaN? | [Limitations and troubleshooting](docs/limitations.md) |
| How do I reproduce an older calculation? | [Historical reproduction](docs/reproducibility.md) |
| Where are APIs, tests and historical reports? | [Developer guide and history](docs/development.md) |
| What FPC functionality exists? | [JET-PLUME appendix](docs/README-JETPLUME.md) |

These pages describe current behavior. Development reports record earlier
milestones and may describe restrictions subsequently removed; they are evidence,
not a substitute for the current manual. The published site may lag this checkout.

<a id="scientific-conventions"></a>

## Scientific conventions

We use Gaussian electromagnetic units and perturbations proportional to
\(\exp(i\mathbf k\cdot\mathbf r-i\omega t)\), with the background magnetic field
along z and \(\mathbf k=(k_\perp,0,k_\parallel)\). Positive
\(\operatorname{Im}\omega\) means growth. Frequencies use the reference gyrofrequency;
input wavevectors and temperatures depend on `normalization='rho_ref'` or
`'d_ref'`. **Do not exchange the two conventions by changing only the flag.**
The [normalization page](docs/normalization.md) supplies the full conversion, including
reference anisotropy, and a verified pair of equivalent inputs.

The default numerical choices are `z_backend='faddeeva'`,
`bessel_mode='adaptive'`, and `heating_convention='corrected'`.
They are independent. [Legacy replay](docs/reproducibility.md) requires explicit
choices and a recorded source revision; fixed harmonics alone are insufficient.

<a id="authors-scientific-references-and-citation"></a>

## Authors, scientific references and citation

Kristopher Klein (kgklein@arizona.edu), Gregory Howes
(gregory-howes@uiowa.edu), and Collin Brown (collin.crbrown@gmail.com).
The kinetic solver descends from Eliot Quataert's solver and Gregory Howes's
Fortran adaptation. Its kinetic equations follow Stix, *Waves in Plasmas*,
Chapter 10, equations 66–73. The cold drifting response follows
[Verscharen & Chandran (2013), Section 2.1](https://arxiv.org/pdf/1212.5192).
The existing selected heating channels are discussed in the
[PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
and [Huang et al. (2024)](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174).

For publications, cite the [PLUME paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2),
acknowledge [the repository](https://github.com/kgklein/PLUME), and record the
version/commit and numerical settings used. The established code-release citation is:

```bibtex
@software{PLUME_2025,
  author = {{Klein}, K. G. and {Howes}, G. G.},
  title = {kgklein/PLUME: Zenodo release},
  month = April,
  year = {2025},
  publisher = {Zenodo},
  version = {v1.0.1},
  doi = {10.5281/zenodo.15215514},
  url = {https://doi.org/10.5281/zenodo.15215514}
}
```

That archived release does not imply that all current features existed in 2025.
The [NASA ADS library](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg)
collects publications using PLUME.

PLUME uses the [BSD 2-Clause license](https://github.com/kgklein/PLUME/blob/main/LICENSE).
The bundled Steven G. Johnson Faddeeva implementation has its own MIT license;
see [third-party provenance](https://github.com/kgklein/PLUME/blob/main/third_party/faddeeva/README.md).
