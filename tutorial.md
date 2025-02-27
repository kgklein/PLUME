# PLUME Tutorial

This is a tutorial for PLUME.
It will guide you through the setting up of some basic input files, the running of the code, and the basic output.
For more details, we refer to the [PLUME Input](input.md) page and the [PLUME Output](output.md) pages.

## Authors

Kristopher Klein   (kgklein@arizona.edu)  
Gregory Howes      (gregory-howes@uiowa.edu)

## Contents

1. Before getting started
2. Installing PLUME
3. Running PLUME


## 1. Before getting started

Before starting with the steps described in this tutorial, we recommend that you familiarise yourself with the code paper.

*code paper here*

You don't need to go through all details, but it is certainly helpful to know what PLUME does and doesn't calculate.
The code paper also explains the numerical techniques used in the code.
We also recommend checking the [Readme](../README.md) file.

## 2. Installing PLUME

This tutorial assumes that you have a working copy of PLUME on your computer, including all the required dependencies. You can find the installation guide [here](../INSTALL.md). Make sure you have a version of PLUME that compiled completely without error messages after typing.

>     make

## 3. Running PLUME

For our first case, consider a simple scan over $k_\parallel \rho_{ref}$ for four solutions, the Slow, Alfven, Fast, and "entropy" (SAFE) modes.

These modes will be identified first by a map scan over a prescribed range of complex frequencies, and then followed along a logirithmic scan of $k_\parallel \rho_{ref}$ values.
```
./plume.e inputs/example/example_map_par.in
```

*Additional cases go here*

