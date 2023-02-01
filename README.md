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
1. Update deltaT for based on CFL for stability
2. Fluids: calculate Darcy term
3. Fluids: PISO Loop (p and U)
4. Platelets Transport: transport via hindered advection-diffusion
5. Biochemicals Transport: transport fluidPhase species via advection-diffusion
6. Platelet Reactions:
  1. Update virtual substance eta
  2. FOR (int m = 0; m < M_rxn; m++ )
    1. React platelets with RK4 solver  
    2. Update platelet boundary conditions 
  3. Update volume fractions for platelets
6. Biochemical Reactions:
  1. FOR (int m = 0; m < M_rxn; m++ )
    1. React biochemcials with RK4 solver  
    2. Update Species fluidPhase boundary conditions 
7. Calculate ADP:
  1. Transport ADP via advection-diffusion
  2. IF sigma_dt has elapsed (e.g. 0.25 s has passed)
    1. Update the source term sigma_release   
8. Write Fields
