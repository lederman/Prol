# PROL
**Code for Computing Generalized Prolate Spheroidal Functions (GPSF)**

The numerical methods for this code are described in: (TBA)

A secondary goal of this project is to develop methodology for cross-validation of key results using both analytical computational validating (here, using Mathematica) and numerical validation.

Code License: GNU General Public License v3.0 (see LICENSE file).

# Project Status
**This code is under development** and it is still being tested. 
The current version does not implement the classic one-dimensional Prolate Spheroidal Wave Functions (PSWF), which we plan to add in the future.
Currently, only a MATLAB implementation is available. We plan to add a more comprehensive implementation in FORTRAN, and interfaces for Julia and Python. 

# Documentation

The numerical methods are described in a paper which is available, along with the associated LaTeX code in the /doc folder of this project. The report is also available on arXiv (TBA).

The code which reproduces the figures in the report is available in the /doc/figures folder

# "Open Source Proof"

To make the analytical expressions more convenient for analytical verification by the readers, we provide experimental Mathematica code that confirms some of the analytical relations which have been obtained in other ways. Where we have not been able to use Mathematica to verify the most general form of the relations that we have derived, we demonstrate some results with specific choice of parameters, which the user can change easily. Furthermore, where possible, we generate the expression in the paper, and some limited pieces of the code, directly from the relations that are verified in the Mathematica code to reduce the possibility of typos and incompatible notation. 
Ultimately, the goal of this experiment is to develop methodology that would allow to confirm key results in papers using both an analytical computational tool and a numerical computational tool, and to verify the compatibility between the expressions. 

# Caveats

* The classic 1-D Prolate Spheroidal Wave Functions (PSWF) have not been implemented in this code yet. 
* The code is still under development, and has not been stress-tested yet. 
* For technical reasons, accuracy testing will only be available with the FORTRAN implementation. 
* The MATLAB implementation relies on the eigenvector decomposition in MATLAB which poses several potential problems in scaling and porting to other languages and other versions of MATLAB. If you attempt to port the code and encounter any loss of precision, please note the comments in that part of the code. This dependency will be removed in future versions.

# Contributing to this Project

We welcome contributions to the code and theoretical background.
The latex source of the report is available in this repository.



