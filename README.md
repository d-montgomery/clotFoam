# clotFoam
## Overview
clotFoam provides a general framework for simulating platelet-mediated coagulation in OpenFOAM.  The solver is based on the platelet aggregation model of Leiderman & Fogelson 2011, with a 12 species coagulation cascade with positive feedback that leads to thrombin generation.  The coagulation model is inspired by Fogelson & Kuharsky 1998.
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
The clotFoam download comes with two tutorials for simulating platelet mediated coagulation.  The rectangle2D case simulates thrombosis in a 2D \[240,60] micron rectangle with an injury length of 90 microns, centered in the middle of the bottom wall of the vessel.  The Hjunction3D case simulates hemostasis in an H-shaped micro fluidic device as described in Schoeman et al.  In both cases, the parameters for the simulation can be edited in the $FOAM_CASE/constan/inputParameters file, and in the $FOAM_CASE/system/controlDict file. To run either of these simulations, navigate back to the main clotFoam directory, then to the desired tutorial directory.  For example tutorials/rectangle2D:

```
$ cd ../tutorials/rectangle2D
Delete any old simulation files (if present):
$ rm -r [1-9]* 0.*
$ blockMesh
$ clotFoam
```

The platelet-mediated coagulation modeled by clotFoam occurs on the real time scale of 10's of minutes.  Therefore, it may take upwards of one day of compute time to simulate clot growth.  

## Parallelization
To run clotFoam in parallel with 6 processors, first edit the decomposeParDict file located in the system directory so that the number of subdomains is 6, and the decomposition method is scotch:
```
numberOfSubdomains 6;

method      scotch;
```
Then following the blockMesh command, decompose the domain and run with 6 processors:
```
$ decomposePar
$ mpirun -np 6 clotFoam -parallel > log &
```

When using an HPC system that utilizes a slurm filesystem, consider using the following outline for the .slurm file for running clotFoam on 2 nodes with a total of 48 cores:
```
#! /bin/bash -x
#SBATCH --job-name="clotFoam_simulation"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
###Multiply nodes and ntasks-per-node and use that for ntasks
#SBATCH --ntasks=48
#SBATCH --time=143:59:00
#SBATCH -o log.slurmFile
#SBATCH -e err.run
# Go to the directory from which our job was launched
cd $SLURM_SUBMIT_DIR

## Get the number of nodes, cores per node, and total number of cores
echo "Number of Nodes = " $SLURM_JOB_NUM_NODES
echo "Number of Cores Per Node = " $SLURM_NTASKS_PER_NODE
echo "Number of Cores = " $SLURM_NTASKS
nodes=$SLURM_JOB_NUM_NODES
cores_per_node=$SLURM_NTASKS_PER_NODE
total_cores=$SLURM_NTASKS

# Create a short JOBID basd on the one provided by the scheduler 
JOBID='echo $SLURM_JOBID'

# Save a copy of our environment and script
cat $0 > script.$JOBID
printenv > env.$JOBID

# Load the necessary modules (these are specific to the HPC system)
module load compilers/gcc/9 
module load mpi/openmpi/gcc/4.1.1
module load apps/openfoam/gcc-openmpi/9
source $OPENFOAM_DIR/etc/bashrc

# Run the job
echo "running job"
srun -N $nodes --ntasks-per-node=$cores_per_node --pty bash
time mpirun -np $total_cores clotFoam -parallel > log 2>&1 
echo "job has finished"               
```

## Algorithm
The solver begins by loading the mesh, reading in constants from constant/inputParameters, reading in fields and boundary conditions from 0/, and initializing the various species objects.  Then the main time-loop is initiated with a dynamically modified time-step based on the maximum Courant number (maxCo) specified in system/controlDict.  First, the solver enters the pressure-velocity loop, where p and U are updated in an iterative sequence known as pressure implicit with splitting of operators (PISO). Next, the platelets and fluid phase biochemicals are transported via advection-diffusion.  Then, the platelets and biochemicals are reacted with one another M times per time step DeltaT. Lastly, the chemical ADP is transported and its source term sigma_release is updated.  The main time-loop iterates until t = endTime, or an error is thrown by the "isSolutionDiverging.H" file.  The algorithm is summarized below:

### clotFoam Algorithm Summary:
* Initialize mesh, constants, fields, and Species objects
* WHILE t < endTime 
  * Update deltaT for based on CFL for stability
  * Fluids: calculate Darcy term
  * Fluids: PISO Loop (p and U)
  * Platelets Transport: 
    * Calcualte hindered velocity flux phi*W:  
      * Interpolate Theta_T to cell faces using downwind scheme
      * Calculate W(Theta_T) for advection
      * Calculate phi*W
    * Calculate hindered diffusion rate Dp*W:
      * Interpolate Theta_T to cell faces using localMax scheme
      * Calculate W(Theta_T) for advection
      * Calculate Dp*W
    * Transport mobile platelets via hindered advection-diffusion
  * Biochemicals Transport: transport all fluidPhase species via advection-diffusion
  * Platelet Reactions:
    * Update virtual substance eta
    * FOR (int m = 0; m < M_rxn; m++ )
      * React platelets with RK4 solver  
    * Update mobile platelet boundary conditions 
    * Update volume fractions for platelets Theta_B and Theta_T
  * Biochemical Reactions:
    * FOR (int m = 0; m < M_rxn; m++ )
       * React biochemcials with RK4 solver  
       * Update Species fluidPhase boundary conditions 
  * Calculate ADP:
    * Transport ADP via advection-diffusion
    * IF sigma_dt has elapsed (e.g. 0.25 s has passed)
      * Update the source term sigma_release   
  * Check if solution is diverging
    * IF Theta_B or Theta_T exceed 1.01: STOP
    * IF pressure p < 0: STOP  
  * Write Fields

## License
OpenFoam, and by extension the clotFoam application, is licensed free and open source under the [GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Acknowledgements
The work was generously supported by grants from the National Science Foundation (NSF DMS-1848221) and the National Institute of Health (NIH R01HL120728 and RO1HL151984). 


## Citing This Work
If you use clotFoam in your work. Please use the following to cite our work:

D. Montgomery, F. Municchi, K. Leiderman, clotFoam: An Open‐Source Framework to Simulate Blood Clot Formation Under Flow, arXiv, 2023, [https://doi.org/10.48550/arXiv.2304.09180](https://doi.org/10.48550/arXiv.2304.09180).


## References
* K. Leiderman and A. L. Fogelson, Grow with the flow: a spatial–temporal model of platelet deposition and blood coagulation under flow. Mathematical Medicine and Biology: a journal of the IMA, 28(1):47–84, 2011. [https://doi.org/10.1093/imammb/dqq005](https://doi.org/10.1093/imammb/dqq005)
* A. L. Fogelson and A. L. Kuharsky. Membrane binding-site density can modulate activation thresholds in enzyme systems. Journal of Theoretical Biology, 193(1):1–18, 1998. [https://doi.org/10.1006/jtbi.1998.0670](https://doi.org/10.1006/jtbi.1998.0670)
* R. M. Schoeman, K. Rana, N. Danes, M. Lehmann, J. A. Di Paola, A. L. Fogelson, K. Leiderman, K.B. Neeves, A microfluidic model of hemostasis sensitive to platelet function and coagulation, Cellular and molecular
bioengineering 10 (2017) 3–15. [https://doi.org/10.1007/s12195-016-0469-0](https://doi.org/10.1007/s12195-016-0469-0)


