---
project: PLUME
project_github: https://github.com/kgklein/PLUME
author: Kristopher Klein, Gregory Howes, Collin Brown
email: kgklein@arizona.edu, gregory-howes@uiowa.edu, collin.crbrown@gmail.com
srd_dir: src
page_dir: docs
output_dir: docs_out
fpp_extensions: f90
base_url: http://www.plume.space
ordered_subpage: INSTALL.md
                 tutorial.md
                 input.md
                 output.md
                 normalization.md
                 MIXED_DISPERSION_REPORT.md
                 FULLY_COLD_DISPERSION_REPORT.md
                 COLD_REFERENCE_MIXED_REPORT.md
                 ADAPTIVE_BESSEL_CACHE_REPORT.md
                 FADDEEVA_BACKEND_REPORT.md
                 MIXED_SCAN_REPORT.md
                 MIXED_EIGENFUNCTION_REPORT.md
                 COLD_FLUID_EIGENFUNCTION_REPORT.md
                 citingpapers.md
                 READ
display: public
         protected
         private
source: false
graph: true
incl_src: false
src_dir:  ./src
favicon: favicon.ico

---

# PLUME: linear waves in kinetic, cold and mixed plasmas

<img src="./page/assets/PLUME_logo.png" alt="PLUME logo" width="200"/>
<img src="./page/assets/qrcode_plume_github.png" alt="Repository QR code" width="200"/>

PLUME is a serial Fortran solver for linear electromagnetic dispersion relations
in homogeneous **nonrelativistic** magnetized plasmas. Populations are drifting
kinetic bi-Maxwellians or collisionless cold fluids, including mixtures with a
cold physical reference species. Parallel equilibrium drifts and kinetic
temperature anisotropy are supported within the documented geometry restrictions.
PLUME does not provide relativistic kinetics, arbitrary distributions or
MPI/OpenMP parallelization.

Start with the [current user manual](./page/index.html),
[installation](./page/INSTALL.html) and [Quick Start](./page/tutorial.html).
The manual is the scientific user entry point; generated procedure lists below
serve developers. New source pages may require staging beyond the historical
workflow; the local verification report records that distinction.

| Workflow or definition | Manual page |
|---|---|
| Every input control | [Input reference](./page/input.html) |
| Physical scales and equivalent inputs | [Normalization](./page/normalization.html) |
| Population models and restrictions | [Species and capabilities](./page/models.html) |
| Roots, frequency maps and scans | [Calculation modes](./page/modes.html) |
| E/B and species density/velocity | [Eigenfunctions](./page/eigenfunctions.html) |
| Signed absorption and complete formal energy | [Heating and wave energy](./page/heating.html) |
| File schemas and Python readers | [Outputs](./page/output.html) |
| Twelve-lesson learning path | [Examples](./page/examples.html) |
| Harmonics, Z and numerical methods | [Numerics](./page/numerics.html) |
| Missing/failed diagnostics | [Limitations](./page/limitations.html) |
| Historical numerical choices | [Reproducibility](./page/reproducibility.html) |
| APIs and historical evidence | [Developer guide](./page/development.html) |

Current outputs include complex mode frequency, electromagnetic polarization,
kinetic/cold density and velocity moments, warm-species signed absorption and
seven selected channels, and corrected wave-energy decomposition. Dispersion
support does not guarantee all supplementary diagnostics: cold kinetic heating
and mixed FPC are unavailable, energy is a qualified weak-mode expression, and
exact-axis restrictions remain. See the capability and limitation pages.

For citation, use the [PLUME code paper](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
and record the source revision/settings. The established release DOI is
[10.5281/zenodo.15215514](https://doi.org/10.5281/zenodo.15215514);
[full citation details](./page/index.html#authors-scientific-references-and-citation)
are in the manual. [Publications using PLUME](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg)
are collected in the existing ADS library.

[JET-PLUME](./page/README-JETPLUME.html) provides historical all-kinetic FPC
reconstructions; no mixed-FPC extension is implied by the new moments/heating.
The optional Python wrapper and readers are described on the output page.

Kristopher Klein, Gregory Howes and Collin Brown maintain PLUME. This manual
update changes documentation only and does not deploy a website.
