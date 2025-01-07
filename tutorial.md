# PLUME Tutorial

This is a tutorial for PLUME.
It will guide you through the setting up of some basic input files, the running of the code, and the basic output.
For more details, we refer to the [PLUME Input](input.md) page, the [PLUME Output](output.md) page, and the [PLUME Documentation](http://plume.space).

## Authors

Kristopher Klein   (kgklein@arizona.edu)  
Gregory Howes      (gregory-howes@uiowa.edu)

## Contents

1. Before getting started
2. Installing PLUME


## 1. Before getting started

Before starting with the steps described in this tutorial, we recommend that you familiarise yourself with the code paper.

[Verscharen, D., Klein, K. G., Chandran, B. D. G., Stevens, M. L., Salem, C. S.,
and Bale, S. D.: ALPS: the Arbitrary Linear Plasma Solver, J. Plasma Phys. 84,
905840403, 2018, doi: 10.1017/S0022377818000739](http://doi.org/10.1017/S0022377818000739)

You don't need to go through all details, but it is certainly helpful to know what ALPS does and doesn't calculate. The code paper also explains the numerical techniques used in the code, and the [ALPS Documentation](http://alps.space) often refers explicitly to equations and sections in the code paper. We also recommend checking the [Readme](../README.md) file.

## 2. Installing PLUME

This tutorial assumes that you have a working copy of ALPS on your computer, including all the required dependencies. You can find the installation guide [here](../INSTALL.md). Make sure you have a version of ALPS that compiled completely without error messages after typing

>     ./configure  
>     make

