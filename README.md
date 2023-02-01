# clotFoam
## Overview
clotFoam provides a general framework for simulating platelet-mediated coagulation in OpenFOAM.  The solver is based on the platelet aggregation model of Leiderman & Fogelson 2011, with an extremely reduced 9 species coagulation cascade that leads to thrombin generation.
The solver is built on the icoFoam code developed by [OpenCFD Ltd.](http://openfoam.com/) to solve the fluids/pressure equations. Target applications for clotFoam include:

* platelet-mediated coagulation
* platelet aggregation
* targeted drug studies
* reactive flows with porous media

## Installation

clotFoam has been developed with the [OpenFoam v9 libraries](https://openfoam.org/version/9/). The code has been tested using the MacOS, Linux, and Ubuntu installations, but should work on any operating system capable of installing OpenFoam. To install the clotFoam solver, first follow the instructions on this page: [OpenFoam v9 Unbuntu Install](https://openfoam.org/download/9-ubuntu/) to install the OpenFoam 9 libraries.  Alternatively, you can find OpenFoam v9 for your desired operating system at the [OpenFoam Download Archive](https://openfoam.org/download/archive/).

After installing OpenFoam v9, navigate to a working folder in a shell terminal, clone the git code repository, and build using OpenFoam v9. <em>Note: MacOS users will need to launch the OpenFoam v9 application using Docker prior to building.</em>

```
$ git clone https://github.com/tomflint22/beamWeldFoam.git clotFoam
$ cd clotFoam/clotFoam
$ wclean
$ wmake
```
