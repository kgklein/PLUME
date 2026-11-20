title: Plasma in a Linear Uniform Magnetized Environment
ordered_subpage: INSTALL.md
                 tutorial.md
                 input.md
                 normalization.md
                 models.md
                 modes.md
                 eigenfunctions.md
                 heating.md
                 output.md
                 examples.md
                 numerics.md
                 limitations.md
                 reproducibility.md
                 development.md
                 README-JETPLUME.md

# User Manual

<img src="./assets/PLUME_logo.png" alt="PLUME logo" width="200"/>
<img src="./assets/qrcode_plume_github.png" alt="Repository QR code" width="200"/>

PLUME solves linear electromagnetic dispersion relations in homogeneous,
nonrelativistic magnetized plasmas, with drifting kinetic bi-Maxwellian and
collisionless cold-fluid populations. This manual documents the physical
conventions, supported calculations, numerical methods, output diagnostics and
limitations of the solver.

For a first calculation, start with [installation](./INSTALL.md) and the
[quick start](./tutorial.md).

## Table of Contents

| Topic | Documentation |
|---|---|
| Building and testing | [Installation](./INSTALL.md) |
| First calculation | [Quick Start](./tutorial.md) |
| Input parameters | [Input reference](./input.md) |
| Physical scales and units | [Normalization](./normalization.md) |
| Kinetic, cold and mixed populations | [Species and models](./models.md) |
| Root finding, maps and scans | [Calculation modes](./modes.md) |
| Fields, density and velocity perturbations | [Eigenfunctions and moments](./eigenfunctions.md) |
| Absorption and formal wave energy | [Heating and wave energy](./heating.md) |
| Output files and Python readers | [Outputs](./output.md) |
| Worked calculations | [Examples](./examples.md) |
| Numerical methods | [Numerics](./numerics.md) |
| Restrictions and troubleshooting | [Limitations](./limitations.md) |
| Historical numerical choices | [Reproducibility](./reproducibility.md) |
| API, tests and development reports | [Developer guide](./development.md) |
| Historical all-kinetic FPC | [JET-PLUME](./README-JETPLUME.md) |

The manual describes current behavior. Development reports document earlier
implementation milestones and may contain restrictions subsequently removed.

<a id="authors-scientific-references-and-citation"></a>

## Authors, scientific references and citation

PLUME is maintained by Kristopher Klein (kgklein@arizona.edu), Gregory Howes
(gregory-howes@uiowa.edu), and Collin Brown (collin.crbrown@gmail.com).
The kinetic solver descends from [Quataert (1998)](https://scixplorer.org/abs/1998ApJ...500..978Q/abstract)
and Gregory Howes's Fortran adaptation. Its kinetic equations follow Stix,
*Waves in Plasmas*, Chapter 10, equations 66–73. The cold drifting response
follows [Verscharen & Chandran (2013), Section 2.1](https://arxiv.org/pdf/1212.5192).
The selected kinetic heating channels are discussed in the
[PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
and [Huang et al. (2024)](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174).

For publications, cite the PLUME code paper, acknowledge the
[repository](https://github.com/kgklein/PLUME), and record the source revision
and numerical settings. The archived code-release citation is:

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

The archived release predates some current capabilities. Publications using
PLUME are collected in the [NASA ADS](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg)
and [NASA SciX](https://scixplorer.org/public-libraries/LR2P1ENiTgO8Z-vBM_35DA)
libraries.

PLUME uses the [BSD 2-Clause license](https://github.com/kgklein/PLUME/blob/main/LICENSE).
The bundled Steven G. Johnson Faddeeva implementation has a separate
[MIT license](https://github.com/kgklein/PLUME/blob/main/third_party/faddeeva/README.md).
