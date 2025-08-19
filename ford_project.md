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

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15215514.svg)](https://doi.org/10.5281/zenodo.15215514)

## Plasma in a Linear Uniform Magnetized Environment

<img src="./page/PLUME_logo.png" alt="drawing" width="200"/>
<img src="./page/qrcode_plume_github.png" alt="drawing" width="200"/>

PLUME is a parallelised numerical code that solves the Vlasov-Maxwell dispersion
relation in hot (even relativistic) magnetised plasma.

If you use the code for a science publication, please provide the code website
[github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements of your publication, and cite the code:

```
@software{PLUME_2025,
  author       = {{Klein}, K. G. and
                  {Howes}, G. G.},
  title        = {kgklein/PLUME: Zenodo release},
  month        = April,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v1.0.1},
  doi          = {10.5281/zenodo.15215514},
  url          = {https://doi.org/10.5281/zenodo.15215514}
}
```
and code paper 
[Klein, K. G., Howes, G. G.,
and Brown, C. R.: PLUME: Plasma in a Linear Uniform Magnetized Environment, RNAAS, 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)

Publications using the PLUME code can be found in our [NASA ADS Library](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg).

For first-time users, we recommend working through our [PLUME Tutorial](./page/tutorial.md).

The key input parameters for PLUME are described on the [PLUME Input](./page/input.md) page.

The output format of PLUME is described on the [PLUME Output](./page/output.md) page.

## Judging Energy Transfer in a - Plasma in a Linear Uniform Magnetized Environment

<img src="./page/Jet-Plume_Logo.svg" alt="JetPlumeLogoDrawing" width="200"/>

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique and linear theory. Please see the [README for JET-PLUME](./page/README-JETPLUME.md).

JET-PLUME is operated in a similar manner to PLUME, either using command line or the python wrapper. See Section 4 above on how to operate and [README for JET-PLUME](./README-JETPLUME.md) for the additional inputs.

Details

If you use this code for a science publication, please provide the same code website as plume
[github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements of your publication and cite the code paper for JET-PLUME: ...

## Python wrapper (linfpclib)

With the creation of JET-PLUME, a wrapper to use PLUME and JET-PLUME in a jupyter notebook was created. It aids in the creation of input files and running of the code with said input files. Please see the [example notebook](./page/examplelinfpc.md) to see the key features of this wrapper. The use of the wrapper is entirely optional. NOTE: Please create and run notebooks from the main directory of the repository. This is necessary because the notebooks import key functions directly from the linfpclib module located there.
