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
display: public
         protected
         private
source: false
graph: true
incl_src: false
src_dir:  ./src
favicon: favicon.ico
preprocess: False

---

[![DOI](https://zenodo.org/badge/243310181.svg)](https://zenodo.org/badge/latestdoi/243310181)
[![build](https://github.com/kgklein/PLUME/actions/workflows/tests.yml/badge.svg)](https://github.com/kgklein/PLUME/actions/workflows/tests.yml)

## Plasma in a Linear Uniform Magnetized Environment

<img src="./page/PLUME_logo.png" alt="drawing" width="200"/>
<img src="./page/qrcode_plume_github.png" alt="drawing" width="200"/>

PLUME is a parallelised numerical code that solves the Vlasov-Maxwell dispersion
relation in hot (even relativistic) magnetised plasma.

If you use the code for a science publication, please provide the code website
[github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements of your publication and cite the code paper:

Publications using the PLUME code can be found in our [NASA ADS Library](https://ui.adsabs.harvard.edu/public-libraries/RWGonkVgRpOaTizvWwsjKg).

---

For first-time users, we recommend working through our [PLUME Tutorial](./page/tutorial.md).

The key input parameters for PLUME are described on the [PLUME Input](./page/input.md) page.

The output format of PLUME is described on the [PLUME Output](./page/output.md) page.

## Judging Energy Transfer in a - Plasma in a Linear Uniform Magnetized Environment

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique and linear theory. 

Details

If you use this code for a science publication, please provide the same code website as plume
[github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements of your publication and cite the code paper for JET-PLUME: ...

## Python wrapper (linfpclib)

With the creation of JET-PLUME, a wrapper to use PLUME and JET-PLUME in a jupyter notebook was created. It aids in the creation of input files and running of the code with said input files. Please see the [example notebook](./page/examplelinfpc.md) to see the key features of this wrapper. The use of the wrapper is entirely optional.