# PLUME: Plasma in a Linear Uniform Magnetized Environment

These are the install instructions for the PLUME code: the Plasma in a Linear Uniform Magnetized Environment.

## Authors

Kristopher Klein   (kgklein@arizona.edu)
Gregory Howes      (gregory-howes@uiowa.edu)

## CONTENTS

1. Requirements and Dependencies
2. Setting up Computer Environments
3. Getting the PLUME Code
4. Installing the PLUME Code
5. Execution of Test Runs

## REQUIREMENTS AND DEPENDENCIES

ALPS has the following requirements:

- A UNIX, Linux, or macOS operating system with a working shell
- GNU make, in some cases it is useful to have the autoconf/automake tools
- Fortran 90 compiler (e.g., gfortran) - we recommend using the latest version
  of the compiler to avoid any surprises in the evaluation.

## SETTING UP COMPUTER ENVIRONMENTS

For Ubuntu and macOS users, the following instructions have proven to be useful.

## GETTING THE PLUME CODE

We recommend pulling the latest version of PLUME from GitHub. For this, go to
the directory where you want to install PLUME and execute the following command:

```
git clone https://github.com/kgklein/PLUME.git
```

Alternatively, you can also go to the website https://https://github.com/kgklein/PLUME

directly and download the source code from there. The advantage of using the git
command is that you can now contribute to the development of the PLUME code. If
you make any changes to the code, GitHub will run automatic tests (via workflows)
to ensure that the changes do not break the code.

## INSTALLING THE PLUME CODE

If all requirements are available, the code can be compiled with the following
commands:

```
make
```

## EXECUTION OF TEST RUNS

PLUME comes with a selection of test runs that cycle through various test
problems. To execute a small set of tests, execute the following shell script:

```
./run_test.sh
```

This script will test the interpolation routine, the routine to generate pre-
described distribution functions, and a simply fast dispersion relation.

## Building the documentation

PLUME uses [Ford](https://forddocs.readthedocs.io/en/latest/) to build its documentation. The documentation is automatically built and deployed to [github.io](https://danielver02.github.io/ALPS/) by the [doc workflow](https://github.com/danielver02/ALPS/blob/master/.github/workflows/doc.yml). To build the documentation locally, follow the [Build documentation](https://github.com/danielver02/ALPS/blob/07a4f8dc996ff76729edeedf5c2a0dc1a5c3028b/.github/workflows/doc.yml#L25-L32) step in the workflow, summarized here:
1. Install `ford` by e.g. `pip install ford`. See [Ford documentation](https://forddocs.readthedocs.io/en/latest/) for details
2. Create a `docs` directory by `mkdir docs`
3. Add a line `title: Readme` to the top of [README.md](./README.md) and copy it to `docs/index.md`
4. Add a line `title: Install` to the top of [INSTALL.md](./INSTALL.md) and copy it to `docs/INSTALL.md`
5. Run `ford ford_project.md`
6. Open `docs/index.html` in a browser to view the documentation

### Adding static pages to documentation

The [README.md](./README.md) and [INSTALL.md](./INSTALL.md) files are added to the [Ford documentation](https://danielver02.github.io/ALPS/) as static pages. You can add more static pages to the documentation by
1. Add the content in a markdown file to the repository.
2. Add a `title: ` line to the beginning of the file and copy it to `docs/` in the [doc workflow](https://github.com/danielver02/ALPS/blob/master/.github/workflows/doc.yml). See steps 3-4 in the previous section, or the [Build documentation](https://github.com/danielver02/ALPS/blob/07a4f8dc996ff76729edeedf5c2a0dc1a5c3028b/.github/workflows/doc.yml#L25-L32) step in the workflow.
3. Add the name of the markdown file as a new line under `ordered_subpage` in [ford_project.md](./ford_project.md)