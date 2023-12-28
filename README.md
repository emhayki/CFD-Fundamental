[![JOSE Article](https://img.shields.io/badge/JOSE-CFD%20Python%3A%20the%2012%20steps%20to%20Navier--Stokes%20equations-<COLOR>.svg)](https://doi.org/10.21105/jose.00021)

This repository presents an implementation of Prof. Lorena A. Barba's ["12 Steps to Navier-Stokes"](https://github.com/barbagroup/CFDPython) tutorial, featuring a methodical approach to understanding and solving the Navier-Stokes equations for fluid flow simulation. Within this repository, you'll find **MATLAB, Python, and C++** code for each of the 12 steps, accompanied by in-depth explanations and references.

## Contents
  - **Steps 1 to 4 in 1D**
   1) Linear Convection
   2) Nonlinear Convection
   3) Diffusion
   4) Burgers' Equation
        
  - **Steps 5 to 10 in 2D**
  5) Linear Convection
  6) Nonlinear Convection
  7) Diffusion
  8) Burgers' Equation
  9) Laplace Equation
  10) Poisson Equation in 2D
   
  - **Steps 11 to 12 - Navier-Stokes Equation in 2D**
  11) Cavity Flow
  12) Channel Flow

 ## C++ and `matplotlibcpp.h`
The C++ codes use [`matplotlibcpp.h`](https://github.com/lava/matplotlib-cpp), a C++ library that provides Matplotlib-like plotting functionality. To use it, include `matplotlibcpp.h` in your C++ project and ensure Python with Matplotlib is installed.

  ## Acknowledgments

- Prof. Lorena A. Barba and the [CFDPython](https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/) project.
- [`matplotlibcpp`](https://github.com/lava/matplotlib-cpp) for C++ plotting capabilities.
