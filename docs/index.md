title: PLUME User Manual
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

# PLUME: linear waves in *cold* and **kinetic** plasmas

<img src="./assets/PLUME_logo.png" alt="PLUME logo" width="200"/>
<img src="./assets/qrcode_plume_github.png" alt="Repository QR code" width="200"/>

PLUME is a FORTRAN dispersion relation solver for linear electromagnetic dispersion relations in homogeneous, magnetized plasmas. 
Populations are treated as drifting hot, anisotropic bi-Maxwellians or cold fluids, with mixtures of hot and cold species allowed. 

PLUME solves linear dispersion relations, outputting
complex mode frequency, 
electromagnetic polarization,
and density and velocity moments. 
For kinetic species, the wave emission or absorption can be calculated, 
including through specified channels, e.g. Landau, transit time, or cyclotron resonances.

New users are encouraged to start with the
[installation](./INSTALL.md) 
and [quick start](./tutorial.md).
This manual serves as the scientific user entry point. 
More detailed API documentation is provided for developers.

| Table of Contents | Manual page |
|---|---|
| Input Parameters | [Input reference](./input.md) |
| Normalization of Physical Scales | [Normalization](./normalization.md) |
| Population models and restrictions | [Species and capabilities](./models.md) |
| Roots, frequency maps and scans | [Calculation modes](./modes.md) |
| E/B and species density/velocity | [Eigenfunctions](./eigenfunctions.md) |
| Signed absorption and complete formal energy | [Heating and wave energy](./heating.md) |
| File schemas and Python readers | [Outputs](./output.md) |
| Twelve-lesson learning path | [Examples](./examples.md) |
| Harmonics, Z and numerical methods | [Numerics](./numerics.md) |
| Missing/failed diagnostics | [Limitations](./limitations.md) |
| Historical numerical choices | [Reproducibility](./reproducibility.md) |
| APIs and historical evidence | [Developer guide](./development.md) |

[JET-PLUME](./README-JETPLUME.md) provides historical all-kinetic FPC
reconstructions. **More Details to follow**


<a id="authors-scientific-references-and-citation"></a>

## Scientific References and Citation

The PLUME kinetic solver descends from Eliot Quataert's solver and Gregory Howes's
FORTRAN adaptation. 
Its kinetic equations follow Stix, *Waves in Plasmas*,
Chapter 10, equations 66–73. 
The selected heating channels are discussed in the
[PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
and [Huang et al. (2024)](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174).
The cold drifting response follows
[Verscharen & Chandran (2013), Section 2.1](https://iopscience.iop.org/article/10.1088/0004-637X/764/1/88/pdf).

For publications, cite the [PLUME paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2),
acknowledge [the repository](https://github.com/kgklein/PLUME), and record the
version/commit and numerical settings used. 
The established code-release citation is:

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

The [NASA SciX library](https://scixplorer.org/public-libraries/LR2P1ENiTgO8Z-vBM_35DA)
collects publications using PLUME. **If we have missed one, please let us know.**

PLUME uses the [BSD 2-Clause license](https://github.com/kgklein/PLUME/blob/main/LICENSE).
The bundled Steven G. Johnson Faddeeva implementation has its own MIT license;
see [third-party provenance](https://github.com/kgklein/PLUME/blob/main/third_party/faddeeva/README.md).
