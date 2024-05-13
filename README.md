# SISC-Paper-2024

[![L2C-CI](https://github.com/spectralDNS/SISC-Paper-2024/actions/workflows/l2c.yml/badge.svg)](https://github.com/spectralDNS/SISC-Paper-2024/actions/workflows/l2c.yml)

Routines for Legendre to Chebyshev (and inverse) transforms. The implemented method is described in the paper 

  * Mikael Mortensen "A faster multipole Legendre-to-Chebyshev transform"

Submitted to SISC. See the preprint [A faster multipole Legendre-to-Chebyshev transform](https://github.com/spectralDNS/SISC-Paper-2024/FMM_paper.pdf).

There are several implementations in the src directory:
  * python - A short, vectorized Python implementation
  * C - An efficient C implementation
  * cython - Wrapping of the fast C code using cython
  * Multiprec - Multiprecision implementations of a direct method using both Python and C++. The multiprecision methods are only used for verification.

In addition there is
  * bin - An executable l2c that can run various tests

# Installation
The code is set up to be compiled with the [meson](https://mesonbuild.com) build system. It should work by cloning this repository and then

    cd src
    meson setup build
    meson install -C build

If you want to install just locally, then use, e.g.,

    meson setup build --prefix=$PWD/build-install

The installation can be tested with

    meson test -C build

# Codespace
Another simple way to test this code is to create a codespace. The environment.yml file in the root folder will then make sure that the codespace creates a conda environment with all necessary dependencies already built. Just press the codespace button and wait awhile for the environment to build. Then enable the environment and run some tests or test the executable `l2c`

     source activate ./venv
     cd src
     meson setup build --prefix=$PWD/build-install --includedir=$CONDA_PREFIX/include --libdir=$CONDA_PREFIX/lib
     meson install -C build
     export PATH=$PWD/build-install/bin:$PATH
     l2c -N1000 -d2 # runs a forward and backward transform and computes the error for an array of length 1000
     meson test -C build # runs all the tests
