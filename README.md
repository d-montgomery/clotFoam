# clotFoam
## Overview
clotFoam provides a general framework for simulating platelet-mediated coagulation in OpenFOAM.  The solver is based on the platelet aggregation model of Leiderman & Fogelson 2011, with an extremely reduced 9 species coagulation cascade that leads to thrombin generation.
The solver is built on the icoFoam code developed by [OpenCFD Ltd.](http://openfoam.com/) to solve the fluids/pressure equations. Target applications for clotFoam include:

* platelet-mediated coagulation
* platelet aggregation
* targeted drug studies
* reactive flows with porous media

## Installation

clotFoam has been developed with the [OpenFoam v9 libraries](https://openfoam.org/version/9/). The code has been tested using the MacOS, Linux, and Ubuntu installations, but should work on any operating system capable of installing OpenFoam. To install the clotFoam solver, first follow the instructions on this page: [OpenFoam v9 Unbuntu Install](https://openfoam.org/download/9-ubuntu/) to install the OpenFoam 9 libraries.  Alternatively, OpenFoam v9 can be downloaded for use with any operating system at the [OpenFoam Download Archive](https://openfoam.org/download/archive/).

After installing OpenFoam v9, navigate to a working folder in a shell terminal, clone the git code repository, and build using OpenFoam v9. <em>Note: MacOS users will need to launch the OpenFoam v9 application using Docker prior to building.</em>

```
$ git clone https://github.com/dmontgomery016/clotFoam.git clotFoam
$ cd clotFoam/clotFoam
$ wclean
$ wmake
```

## Tutorial cases
The clotFoam download comes with a 2D base case for simulating platelet mediated coagulation in a \[240,60] micron rectangle. The injury length is set to 90 microns, centered in the middle of the bottom wall of the vessel.  To run this simulation, navigate back to the main clotFoam directory, then to the baseCaseClotFoam_Rectangle directory:

```
$ cd ../baseCaseClotFoam_Rectangle
Delete any old simulation files (if present):
$ rm -r [1-9]* 0.*
$ blockMesh
$ clotFoam
```

The platelet-mediated coagulation modeled by clotFoam occurs on the real time scale of 10s of minutes.  Therefore, it may take upwards of one day of compute time to simulate clot growth.  

## Parallelization
To run clotFoam in parallel with P processors, first edit the decomposeParDict file located in the system directory so that the number of subdomains is P, and the decomposition method is scotch:
```
numberOfSubdomains P;

method      scotch;
```
Then following the blockMesh command, decompose the domain and run with P processors:
```
$ decomposePar
$ mpirun -np P clotFoam -parallel > log &
```

