---
project: PLUME
project_github: https://github.com/kgklein/PLUME
author: Kristopher Klein, Gregory Howes
email: kgklein@arizona.edu, gregory-howes@uiowa.edu
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
[github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements of your publication and cite the code paper:

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

Publications using the PLUME code can be found in our [NASA ADS Library](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg).

---

For first-time users, we recommend working through our [PLUME Tutorial](./page/tutorial.md).

The key input parameters for PLUME are described on the [PLUME Input](./page/input.md) page.

The output format of PLUME is described on the [PLUME Output](./page/output.md) page.
